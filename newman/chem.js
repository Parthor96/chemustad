/*
 * newman/chem.js  —  window.NNChem (Newman Navigator, chemistry input layer)
 *
 * Turns what a student types (a name, a molecular formula, a condensed formula,
 * or SMILES) or draws (via NNSketch -> fromGraph) into a heavy-atom molecular
 * graph with implicit hydrogens, rings, and substituent "group" labels.
 * Contract: newman/SPEC.md sections 2 and 3. No dependencies and no random numbers.
 *
 * Deviations / resolved ambiguities (signatures are unchanged):
 *  - EL.S.valence is [2, 4, 6] and EL.P.valence is [3, 5] (the spec lists [2] / [3])
 *    so SMILES like CS(=O)C get standard implicit-H counts. The first entry is
 *    still the "normal" valence used by the condensed parser and the sketch.
 *  - Every successful parse gets a suggested defaultBond, not only presets: the
 *    rotatable heavy-atom single bond with the most heavy atoms on its smaller
 *    side (ties: C–C first, then lowest indices). Falls back to a ring bond, else null.
 *  - SMILES typed exactly as a preset / known-name SMILES (e.g. "CCCC") gets that
 *    name, presetId and defaultBond; source stays 'smiles'.
 *  - Besides the 28 presets, a small dictionary of extra names (methanol,
 *    2-propanol, diethyl ether, toluene, ...) is accepted by name. The typed name
 *    "2,3-dibromobutane" returns choices (meso or 2R,3R).
 *  - Explicit [H] atoms (SMILES) and 'H' vertices (fromGraph) are folded into
 *    their heavy atom's hCount so heavy atoms always occupy 0..nHeavy-1; a chiral
 *    neighbour slot that pointed at such an H becomes -1.
 *  - Error results may carry an extra `atom` field (heavy index) for valence
 *    errors, so the sketch can ring the offending vertex.
 *  - Condensed-formula ambiguities (is a bare O in "CH3COCH3" a chain O or a C=O?)
 *    are resolved by trying chain-first, then C=O, keeping the first assignment
 *    where every atom ends with a legal valence.
 *  - A formula-shaped input that is not in the isomer table (CH3Cl, CH2Cl2) is
 *    retried as a condensed formula before giving up.
 *  - The CnH2n+2 fallback (n = 7..10) names its isomers (heptane, 2-methylhexane, ...).
 */
(function (root) {
  'use strict';

  /* ------------------------------------------------------------------ */
  /* Element table                                                        */
  /* ------------------------------------------------------------------ */
  const EL = {
    H:  { color3d: 0xf5f0e6, color2d: '#fbf3e2', r: 0.30, vdw: 1.20, valence: [1] },
    C:  { color3d: 0x3a3a3a, color2d: '#ffb300', r: 0.48, vdw: 1.70, valence: [4] },
    N:  { color3d: 0x3b6fe0, color2d: '#6f9bff', r: 0.46, vdw: 1.55, valence: [3] },
    O:  { color3d: 0xe53935, color2d: '#ff6b61', r: 0.46, vdw: 1.52, valence: [2] },
    F:  { color3d: 0x9be15d, color2d: '#9be15d', r: 0.40, vdw: 1.47, valence: [1] },
    Cl: { color3d: 0x2ec27e, color2d: '#2ec27e', r: 0.50, vdw: 1.75, valence: [1] },
    Br: { color3d: 0xa0412d, color2d: '#e0875f', r: 0.53, vdw: 1.85, valence: [1] },
    I:  { color3d: 0x7e3fbf, color2d: '#b784f0', r: 0.58, vdw: 1.98, valence: [1] },
    S:  { color3d: 0xf2c200, color2d: '#f2c200', r: 0.54, vdw: 1.80, valence: [2, 4, 6] },
    B:  { color3d: 0xffa07a, color2d: '#ffa07a', r: 0.44, vdw: 1.92, valence: [3] },
    P:  { color3d: 0xff8c00, color2d: '#ff8c00', r: 0.54, vdw: 1.80, valence: [3, 5] }
  };
  const EL_NAME = {
    H: 'Hydrogen', C: 'Carbon', N: 'Nitrogen', O: 'Oxygen', F: 'Fluorine', Cl: 'Chlorine',
    Br: 'Bromine', I: 'Iodine', S: 'Sulfur', B: 'Boron', P: 'Phosphorus'
  };
  const HALOGEN = { F: true, Cl: true, Br: true, I: true };
  const HEAVY_MAX = 30;

  /* Allowed total valences (bond orders + H) for an element with a charge. */
  const CHARGED_VALENCE = {
    'N+1': [4], 'O+1': [3], 'O-1': [1], 'N-1': [2], 'C+1': [3], 'C-1': [3],
    'S+1': [3], 'S-1': [1], 'B-1': [4], 'P+1': [4]
  };
  function valenceList(el, charge) {
    const e = EL[el];
    if (!e) return [];
    charge = charge | 0;
    if (!charge) return e.valence;
    const key = el + (charge > 0 ? '+' : '-') + Math.abs(charge);
    if (CHARGED_VALENCE[key]) return CHARGED_VALENCE[key];
    if (HALOGEN[el] && charge === -1) return [0];
    return [Math.max(0, e.valence[0] - Math.abs(charge))];
  }
  /* Smallest allowed valence >= used; -1 if used exceeds every allowed valence. */
  function fitValence(el, charge, used) {
    const list = valenceList(el, charge);
    for (let k = 0; k < list.length; k++) if (list[k] >= used) return list[k];
    return -1;
  }
  function maxValence(el, charge) {
    const list = valenceList(el, charge);
    return list.length ? list[list.length - 1] : 0;
  }

  /* ------------------------------------------------------------------ */
  /* Small graph helpers                                                  */
  /* ------------------------------------------------------------------ */
  function neighbors(mol, i) {
    const out = [];
    for (let k = 0; k < mol.bonds.length; k++) {
      const b = mol.bonds[k];
      if (b.a === i) out.push(b.b); else if (b.b === i) out.push(b.a);
    }
    return out.sort((x, y) => x - y);
  }
  function bondBetween(mol, i, j) {
    for (let k = 0; k < mol.bonds.length; k++) {
      const b = mol.bonds[k];
      if ((b.a === i && b.b === j) || (b.a === j && b.b === i)) return k;
    }
    return -1;
  }
  function heavyNeighbors(mol, i) {
    return neighbors(mol, i).filter(n => mol.atoms[n].el !== 'H');
  }
  function orderSums(nAtoms, bonds) {
    const s = new Array(nAtoms).fill(0);
    for (const b of bonds) { s[b.a] += b.order; s[b.b] += b.order; }
    return s;
  }
  /* Number of connected pieces among atoms 0..n-1. */
  function countPieces(n, bonds) {
    if (!n) return 0;
    const parent = []; for (let i = 0; i < n; i++) parent.push(i);
    const find = x => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    for (const b of bonds) { const ra = find(b.a), rb = find(b.b); if (ra !== rb) parent[ra] = rb; }
    let c = 0; for (let i = 0; i < n; i++) if (find(i) === i) c++;
    return c;
  }
  function hybrid(mol, i) {
    let dbl = 0, tri = 0;
    for (const b of mol.bonds) {
      if (b.a !== i && b.b !== i) continue;
      if (b.order === 2) dbl++; else if (b.order === 3) tri++;
    }
    if (tri || dbl >= 2) return 'sp';
    if (dbl) return 'sp2';
    return 'sp3';
  }
  /* Display numbering (integration fix). Students expect carbons numbered along the main
   * chain (1,2-dichloroethane's C–C bond is C1–C2, meso-2,3-dibromobutane's centre is
   * C2–C3), not by position in the input string. Carbons: the main chain first (largest
   * carbon ring, else the longest carbon chain), direction and start chosen for the most
   * then lowest substituent locants, remaining carbons in BFS order. Other elements are
   * numbered per element (Cl1, Cl2), and a heteroatom that occurs once has no number (O). */
  const numCache = typeof WeakMap !== 'undefined' ? new WeakMap() : null;
  function displayNumbers(mol) {
    if (numCache && numCache.has(mol)) {
      const c = numCache.get(mol);
      if (c.n === mol.atoms.length && c.nb === mol.bonds.length) return c.num;
    }
    const n = mol.atoms.length;
    const heavy = [];
    for (let i = 0; i < n; i++) if (mol.atoms[i].el !== 'H') heavy.push(i);
    const adj = mol.atoms.map(() => []);
    const multi = new Array(n).fill(0);
    for (const b of mol.bonds) {
      if (mol.atoms[b.a].el === 'H' || mol.atoms[b.b].el === 'H') continue;
      adj[b.a].push(b.b); adj[b.b].push(b.a);
      if (b.order > 1) { multi[b.a]++; multi[b.b]++; }
    }
    const isC = (i) => mol.atoms[i].el === 'C';
    // locant vector for an ordered main chain: positions of substituents and multiple bonds
    const score = (path) => {
      const set = new Set(path), loc = [];
      path.forEach((a, k) => {
        adj[a].forEach((x) => { if (!set.has(x)) loc.push(k + 1); });
        if (multi[a]) loc.push(k + 1); // unsaturation counts at its atom
      });
      return loc.sort((x, y) => x - y);
    };
    const better = (la, lb) => { // true if locant list la beats lb
      if (la.length !== lb.length) return la.length > lb.length;
      for (let k = 0; k < la.length; k++) if (la[k] !== lb[k]) return la[k] < lb[k];
      return false;
    };
    let main = null, mainScore = null;
    const consider = (path) => {
      const sc = score(path);
      if (!main || path.length > main.length || (path.length === main.length && better(sc, mainScore))) { main = path; mainScore = sc; }
    };
    // largest carbon-containing ring (carbons only along it)
    const rings = (mol.rings || []).slice().sort((x, y) => y.filter(isC).length - x.filter(isC).length);
    if (rings.length && rings[0].filter(isC).length >= 3) {
      const best = rings[0].filter(isC).length;
      rings.filter((q) => q.filter(isC).length === best).forEach((rr) => {
        const L = rr.length;
        for (let s0 = 0; s0 < L; s0++) for (const dir of [1, -1]) {
          const path = [];
          for (let k = 0; k < L; k++) path.push(rr[((s0 + dir * k) % L + L) % L]);
          const cp = path.filter(isC);
          if (main && main.length > cp.length) continue;
          // locants over carbon positions only
          const set = new Set(path), loc = [];
          cp.forEach((a, k) => { adj[a].forEach((x) => { if (!set.has(x)) loc.push(k + 1); }); });
          loc.sort((x, y) => x - y);
          if (!main || cp.length > main.length || (cp.length === main.length && better(loc, mainScore))) { main = cp; mainScore = loc; }
        }
      });
    } else {
      // longest simple path through carbons (carbon subgraph is a forest when acyclic)
      const cs = heavy.filter(isC);
      const cadj = (i) => adj[i].filter(isC);
      for (const s0 of cs) {
        if (cadj(s0).length > 1) continue; // chain ends only
        const prev = new Map([[s0, -1]]), q = [s0], order = [];
        while (q.length) { const x = q.shift(); order.push(x); for (const y of cadj(x)) if (!prev.has(y)) { prev.set(y, x); q.push(y); } }
        for (const e of order) {
          if (e === s0 && cs.length > 1) continue;
          const path = []; let c = e;
          while (c !== -1) { path.push(c); c = prev.get(c); }
          path.reverse();
          consider(path);
        }
      }
      if (!main && cs.length) main = [cs[0]];
    }
    const num = new Array(n).fill(0);
    let next = 1;
    (main || []).forEach((a) => { num[a] = next++; });
    // remaining carbons: BFS outward from the main chain, lower input index first
    const seen = new Set(main || []);
    let frontier = (main || []).slice();
    while (frontier.length) {
      const nf = [];
      frontier.forEach((x) => adj[x].slice().sort((a, b) => a - b).forEach((y) => {
        if (seen.has(y)) return; seen.add(y); nf.push(y);
        if (isC(y)) num[y] = next++;
      }));
      frontier = nf;
    }
    heavy.forEach((i) => { if (isC(i) && !num[i]) num[i] = next++; });
    // other elements: per-element counters in main-chain-first order
    const byEl = {};
    const orderAll = heavy.slice().sort((a, b) => {
      const da = dist0(a), db = dist0(b); return da - db || a - b;
    });
    function dist0(i) { // first carbon neighbour's number, to follow the chain
      const cn = adj[i].filter(isC).map((x) => num[x]);
      return cn.length ? Math.min(...cn) : 1e6;
    }
    const counts = {};
    heavy.forEach((i) => { if (!isC(i)) counts[mol.atoms[i].el] = (counts[mol.atoms[i].el] || 0) + 1; });
    orderAll.forEach((i) => { if (!isC(i)) { const e = mol.atoms[i].el; byEl[e] = (byEl[e] || 0) + 1; num[i] = counts[e] > 1 ? byEl[e] : 0; } });
    if (numCache) numCache.set(mol, { n, nb: mol.bonds.length, num });
    return num;
  }
  function atomLabel(mol, i) {
    const a = mol.atoms[i];
    if (!a) return '?';
    if (a.el === 'H') return 'H';
    let k;
    try { k = displayNumbers(mol)[i]; } catch (e) { k = i + 1; }
    return a.el + (k ? k : (a.el === 'C' ? i + 1 : ''));
  }
  function atomNumber(mol, i) { try { return displayNumbers(mol)[i] || 0; } catch (e) { return i + 1; } }
  function bondLabel(mol, a, b) {
    // carbons first, then by display number (C1–C2, C2–O)
    const key = (i) => (mol.atoms[i].el === 'C' ? 0 : 1000) + atomNumber(mol, i);
    const lo = key(a) <= key(b) ? a : b, hi = lo === a ? b : a;
    return atomLabel(mol, lo) + '–' + atomLabel(mol, hi);
  }

  /* ------------------------------------------------------------------ */
  /* Ring perception (heavy atoms): shortest cycle per ring bond, then a  */
  /* GF(2) independence filter down to E - V + C rings (approximate SSSR) */
  /* ------------------------------------------------------------------ */
  function perceiveRings(atoms, bonds) {
    const heavy = i => atoms[i].el !== 'H';
    const adj = atoms.map(() => []);
    let nE = 0, nV = 0;
    bonds.forEach((b, e) => {
      b.inRing = false;
      if (heavy(b.a) && heavy(b.b)) {
        adj[b.a].push({ n: b.b, e }); adj[b.b].push({ n: b.a, e }); nE++;
      }
    });
    adj.forEach(l => l.sort((x, y) => x.n - y.n));
    const heavyIdx = [];
    atoms.forEach((a, i) => { if (heavy(i)) { heavyIdx.push(i); nV++; } });
    const heavyBonds = bonds.filter(b => heavy(b.a) && heavy(b.b));
    // components among heavy atoms (heavy atoms are 0..nV-1 by contract, but be safe)
    const remap = new Map(); heavyIdx.forEach((h, k) => remap.set(h, k));
    const comps = countPieces(nV, heavyBonds.map(b => ({ a: remap.get(b.a), b: remap.get(b.b) })));
    const nRings = nE - nV + comps;
    if (nRings <= 0) return [];

    const cands = []; const seen = new Set();
    bonds.forEach((bd, e) => {
      if (!heavy(bd.a) || !heavy(bd.b)) return;
      // BFS from a to b without using edge e
      const prev = new Map(); prev.set(bd.a, null);
      const q = [bd.a]; let found = false;
      while (q.length && !found) {
        const x = q.shift();
        for (const { n, e: e2 } of adj[x]) {
          if (e2 === e || prev.has(n)) continue;
          prev.set(n, { from: x, e: e2 });
          if (n === bd.b) { found = true; break; }
          q.push(n);
        }
      }
      if (!found) return;
      bd.inRing = true;
      const path = [bd.b], edges = [e];
      let cur = bd.b;
      while (prev.get(cur)) { const p = prev.get(cur); edges.push(p.e); path.push(p.from); cur = p.from; }
      const key = edges.slice().sort((x, y) => x - y).join(',');
      if (seen.has(key)) return;
      seen.add(key);
      cands.push({ path: path.reverse(), edges, key });
    });
    cands.sort((x, y) => x.edges.length - y.edges.length || (x.key < y.key ? -1 : x.key > y.key ? 1 : 0));

    const basis = []; const rings = [];
    for (const c of cands) {
      if (rings.length >= nRings) break;
      const v = new Uint8Array(bonds.length);
      c.edges.forEach(e => { v[e] = 1; });
      for (const b of basis) if (v[b.pivot]) for (let k = 0; k < v.length; k++) v[k] ^= b.vec[k];
      const pivot = v.indexOf(1);
      if (pivot < 0) continue;
      basis.push({ vec: v, pivot });
      rings.push(normalizeCycle(c.path));
    }
    return rings;
  }
  /* Rotate a cycle so the smallest index is first and the walk heads to its smaller neighbour. */
  function normalizeCycle(path) {
    const n = path.length;
    let m = 0; for (let k = 1; k < n; k++) if (path[k] < path[m]) m = k;
    let out = path.slice(m).concat(path.slice(0, m));
    if (n > 2 && out[n - 1] < out[1]) out = [out[0]].concat(out.slice(1).reverse());
    return out;
  }

  /* ------------------------------------------------------------------ */
  /* perceive / addHydrogens / formula                                    */
  /* ------------------------------------------------------------------ */
  function perceive(mol) {
    const atoms = mol.atoms;
    let nHeavy = 0;
    for (const a of atoms) if (a.el !== 'H') nHeavy++;
    mol.nHeavy = nHeavy;
    const sums = orderSums(atoms.length, mol.bonds);
    const hNbr = new Array(atoms.length).fill(0);
    for (const b of mol.bonds) {
      if (atoms[b.a].el === 'H') hNbr[b.b]++;
      if (atoms[b.b].el === 'H') hNbr[b.a]++;
    }
    atoms.forEach((a, i) => {
      if (a.charge == null) a.charge = 0;
      if (a.el === 'H') { a.hCount = 0; a.hTotal = 0; return; }
      if (a.hCount == null || !(a.hCount >= 0)) {
        const v = fitValence(a.el, a.charge, sums[i]);
        a.hCount = v < 0 ? 0 : v - sums[i];
      }
      a.hTotal = a.hCount + hNbr[i];
    });
    mol.rings = perceiveRings(atoms, mol.bonds);
    if (mol.explicitH == null) mol.explicitH = atoms.length > nHeavy;
    return mol;
  }

  function cloneMol(mol) {
    return {
      atoms: mol.atoms.map(a => {
        const c = Object.assign({}, a);
        if (a.chiralNbrs) c.chiralNbrs = a.chiralNbrs.slice();
        return c;
      }),
      bonds: mol.bonds.map(b => (b.stereo ? { a: b.a, b: b.b, order: b.order, inRing: !!b.inRing, stereo: Object.assign({}, b.stereo) }
        : { a: b.a, b: b.b, order: b.order, inRing: !!b.inRing })),
      nHeavy: mol.nHeavy,
      explicitH: !!mol.explicitH,
      rings: (mol.rings || []).map(r => r.slice()),
      name: mol.name == null ? null : mol.name,
      source: mol.source,
      input: mol.input == null ? '' : mol.input,
      defaultBond: mol.defaultBond ? mol.defaultBond.slice() : null,
      presetId: mol.presetId == null ? null : mol.presetId
    };
  }

  function addHydrogens(mol) {
    const m = cloneMol(mol);
    if (mol.explicitH) return m;
    const nH = m.atoms.filter(a => a.el !== 'H').length;
    for (let i = 0; i < nH; i++) {
      const a = m.atoms[i];
      const k = a.hCount | 0;
      let first = -1;
      for (let j = 0; j < k; j++) {
        const idx = m.atoms.length;
        if (first < 0) first = idx;
        m.atoms.push({ el: 'H', charge: 0, hCount: 0, hTotal: 0, parent: i });
        m.bonds.push({ a: i, b: idx, order: 1, inRing: false });
      }
      a.hTotal = (a.hTotal == null ? 0 : a.hTotal);
      if (a.hTotal < k) a.hTotal = k;
      a.hCount = 0;
      if (a.chiralNbrs && first >= 0) {
        const p = a.chiralNbrs.indexOf(-1);
        if (p >= 0) a.chiralNbrs[p] = first;
      }
    }
    m.explicitH = true;
    m.nHeavy = nH;
    return m;
  }

  function elementCounts(mol) {
    const c = {};
    for (const a of mol.atoms) {
      c[a.el] = (c[a.el] || 0) + 1;
      if (a.el !== 'H' && a.hCount) c.H = (c.H || 0) + a.hCount;
    }
    return c;
  }
  function hill(counts) {
    const keys = Object.keys(counts).filter(k => counts[k] > 0);
    const part = k => k + (counts[k] > 1 ? counts[k] : '');
    let out = '';
    if (counts.C) {
      out += part('C');
      if (counts.H) out += part('H');
      keys.filter(k => k !== 'C' && k !== 'H').sort().forEach(k => { out += part(k); });
    } else {
      keys.sort().forEach(k => { out += part(k); });
    }
    return out;
  }
  function formula(mol) {
    let q = 0;
    for (const a of mol.atoms) q += a.charge | 0;
    return hill(elementCounts(mol)) + chargeText(q);
  }

  /* ------------------------------------------------------------------ */
  /* Group labels                                                         */
  /* ------------------------------------------------------------------ */
  const HETERO_SCORE = { I: 2.9, Br: 2.8, S: 2.7, Cl: 2.6, N: 2.3, O: 2.2, F: 2.0 };

  function sharesRing(mol, i, j) {
    const rings = mol.rings || [];
    for (const r of rings) if (r.indexOf(i) >= 0 && r.indexOf(j) >= 0) return true;
    return false;
  }
  function inAnyRing(mol, i) {
    const rings = mol.rings || [];
    for (const r of rings) if (r.indexOf(i) >= 0) return true;
    return false;
  }
  function hPart(n) { return n <= 0 ? '' : n === 1 ? 'H' : 'H' + n; }
  function chargeText(q) {
    if (!q) return '';
    const s = q > 0 ? '+' : '−';
    return Math.abs(q) > 1 ? Math.abs(q) + s : s;
  }
  function hTotalOf(mol, i) {
    const a = mol.atoms[i];
    if (a.hTotal != null) return a.hTotal;
    let n = a.hCount | 0;
    for (const j of neighbors(mol, i)) if (mol.atoms[j].el === 'H') n++;
    return n;
  }
  /* Heavy atoms reachable from s without crossing `from` (s included). */
  function subtree(mol, s, from) {
    const seen = new Set([s]); const q = [s];
    while (q.length) {
      const x = q.shift();
      for (const n of heavyNeighbors(mol, x)) {
        if (n === from || seen.has(n)) continue;
        seen.add(n); q.push(n);
      }
    }
    return Array.from(seen);
  }

  /* Label for a substituent that starts at ring atom s (s's ring does not contain `from`). */
  function ringSubLabel(mol, s) {
    const rings = (mol.rings || []).filter(r => r.indexOf(s) >= 0).sort((x, y) => x.length - y.length);
    const r = rings[0];
    const a = mol.atoms[s];
    const fallback = a.el + hPart(hTotalOf(mol, s)) + 'R';
    if (!r || rings.length > 1) return fallback;
    const allC = r.every(i => mol.atoms[i].el === 'C');
    // substituents on ring atoms other than the attachment bond
    let extra = 0;
    for (const i of r) {
      for (const n of heavyNeighbors(mol, i)) if (r.indexOf(n) < 0) extra++;
    }
    extra -= 1; // the bond back toward `from`
    let dbl = 0;
    for (let k = 0; k < r.length; k++) {
      const bi = bondBetween(mol, r[k], r[(k + 1) % r.length]);
      if (bi >= 0 && mol.bonds[bi].order === 2) dbl++;
    }
    if (allC && r.length === 6 && dbl === 3) return extra === 0 ? 'Ph' : 'Ar';
    if (allC && dbl === 0 && extra === 0) return 'C' + r.length + 'H' + (2 * r.length - 1);
    return fallback;
  }

  /* Condensed label for the group starting at s, walking away from `from`. */
  function condLabel(mol, s, from, visited) {
    visited.add(s);
    const a = mol.atoms[s];
    if (a.el === 'H') return 'H';
    if (from != null && inAnyRing(mol, s) && !sharesRing(mol, s, from)) return ringSubLabel(mol, s);
    const chg = chargeText(a.charge);
    const H = hPart(hTotalOf(mol, s));
    const kids = heavyNeighbors(mol, s).filter(n => n !== from && !visited.has(n)).map(n => {
      const bi = bondBetween(mol, s, n);
      return { idx: n, el: mol.atoms[n].el, order: bi >= 0 ? mol.bonds[bi].order : 1 };
    });
    if (HALOGEN[a.el] && !kids.length) return a.el + chg;
    if (!kids.length) return a.el + H + chg;
    kids.forEach(k => { k.label = condLabel(mol, k.idx, s, visited); });

    let out = a.el + H;
    // carbonyl-type =O / =S written inline right after the atom: CHO, COOH, COCH3
    let oxo = kids.filter(k => k.order === 2 && (k.el === 'O' || k.el === 'S') && k.label === k.el);
    let rest = kids.filter(k => oxo.indexOf(k) < 0);
    if (oxo.length === 1) {
      const est = rest.find(k => k.order === 1 && k.el === 'O' && k.label.length > 1 &&
        k.label !== 'OH' && k.label.charAt(1) !== '−');
      if (est) { // ester: CO2CH3
        out += 'O2' + est.label.slice(1);
        rest = rest.filter(k => k !== est);
        oxo = [];
      }
    }
    oxo.forEach(k => { out += k.label; });
    if (!rest.length) return out + chg;

    const sym = o => (o === 2 ? '=' : o === 3 ? '≡' : '');
    const allSame = rest.length > 1 && rest.every(k => k.label === rest[0].label && k.order === rest[0].order);
    let tail = null;
    if (!allSame) {
      const multi = rest.filter(k => k.order > 1);
      if (multi.length) tail = multi[multi.length - 1];
      else {
        const sorted = rest.slice().sort((x, y) => x.label.length - y.label.length ||
          (x.label < y.label ? -1 : x.label > y.label ? 1 : 0));
        tail = sorted[sorted.length - 1];
      }
    }
    const groups = [];
    for (const k of rest) {
      if (k === tail) continue;
      const key = sym(k.order) + k.label;
      const g = groups.find(x => x.key === key);
      if (g) g.n++; else groups.push({ key, label: k.label, order: k.order, n: 1 });
    }
    groups.sort((x, y) => x.key.length - y.key.length || (x.key < y.key ? -1 : x.key > y.key ? 1 : 0));
    for (const g of groups) {
      const cnt = g.n > 1 ? String(g.n) : '';
      if (g.order === 1 && HALOGEN[g.label]) out += g.label + cnt;
      else out += '(' + g.key + ')' + cnt;
    }
    if (tail) out += sym(tail.order) + tail.label;
    return out + chg;
  }

  /* Short names for long alkyl labels (Pr, Bu, sBu, iBu), else first atom + R. */
  function shortLabel(mol, s, from) {
    const sub = subtree(mol, s, from);
    const alkyl = sub.every(i => mol.atoms[i].el === 'C' && hybrid(mol, i) === 'sp3' && !inAnyRing(mol, i));
    if (alkyl) {
      const kids = x => heavyNeighbors(mol, x).filter(n => n !== from && sub.indexOf(n) >= 0 && n !== x);
      const k0 = heavyNeighbors(mol, s).filter(n => n !== from);
      if (sub.length === 2) return 'Et';
      if (sub.length === 3) return k0.length === 2 ? 'iPr' : 'Pr';
      if (sub.length === 4) {
        if (k0.length === 3) return 'tBu';
        if (k0.length === 2) return 'sBu';
        const c1 = k0[0];
        const k1 = heavyNeighbors(mol, c1).filter(n => n !== s);
        return k1.length === 2 ? 'iBu' : 'Bu';
      }
      void kids;
    }
    const a = mol.atoms[s];
    return a.el + hPart(hTotalOf(mol, s)) + 'R';
  }

  function groupInfo(mol, s, from) {
    const a = mol.atoms[s];
    if (!a) return { label: '?', cls: 'H', score: 0, el: '?', ring: false, hasH: false };
    if (a.el === 'H') return { label: 'H', cls: 'H', score: 0, el: 'H', ring: false, hasH: false };
    const ring = from != null && sharesRing(mol, s, from);
    const hTot = hTotalOf(mol, s);
    const kids = heavyNeighbors(mol, s).filter(n => n !== from);
    const hyb = hybrid(mol, s);
    let cls, score;
    if (a.el === 'C' || a.el === 'B' || a.el === 'P') {
      if (a.el !== 'C' || hyb !== 'sp3' || kids.length <= 1) cls = 'Me';
      else cls = kids.length === 2 ? 'iPr' : 'tBu';
    } else {
      cls = a.el;
    }
    if (cls === 'tBu') score = 5;
    else if (cls === 'iPr') score = 4;
    else if (cls === 'Me') score = (a.el === 'C' && hyb !== 'sp3') ? 3.05 : 3 + 0.01 * subtree(mol, s, from).length;
    else score = HETERO_SCORE[cls] != null ? HETERO_SCORE[cls] : 2.0;

    let label;
    if (ring) {
      label = a.el + hPart(hTot) + chargeText(a.charge);
    } else {
      label = condLabel(mol, s, from, new Set(from != null ? [from] : []));
      if (label.length > 8) label = shortLabel(mol, s, from);
    }
    if (from != null) {
      const bi = bondBetween(mol, s, from);
      const o = bi >= 0 ? mol.bonds[bi].order : 1;
      if (o === 2) label = '=' + label; else if (o === 3) label = '≡' + label;
    }
    const hasH = (a.el === 'O' || a.el === 'N' || a.el === 'S') && hTot > 0;
    return { label, cls, score, el: a.el, ring, hasH };
  }
  function groupLabel(mol, s, from) { return groupInfo(mol, s, from).label; }

  /* ------------------------------------------------------------------ */
  /* Presets and names                                                    */
  /* ------------------------------------------------------------------ */
  function P(id, name, aliases, smiles, defaultBond, group) {
    return { id, name, aliases, smiles, defaultBond, group };
  }
  const presets = [
    P('ethane', 'ethane', [], 'CC', [0, 1], 'Alkanes'),
    P('propane', 'propane', [], 'CCC', [0, 1], 'Alkanes'),
    P('butane', 'butane', ['n-butane'], 'CCCC', [1, 2], 'Alkanes'),
    P('isobutane', '2-methylpropane', ['isobutane'], 'CC(C)C', [0, 1], 'Alkanes'),
    P('pentane', 'pentane', ['n-pentane'], 'CCCCC', [1, 2], 'Alkanes'),
    P('2-methylbutane', '2-methylbutane', ['isopentane'], 'CC(C)CC', [1, 3], 'Alkanes'),
    P('2,2-dimethylbutane', '2,2-dimethylbutane', [], 'CC(C)(C)CC', [1, 4], 'Alkanes'),
    P('2,3-dimethylbutane', '2,3-dimethylbutane', [], 'CC(C)C(C)C', [1, 3], 'Alkanes'),
    P('hexane', 'hexane', ['n-hexane'], 'CCCCCC', [2, 3], 'Alkanes'),
    P('chloroethane', 'chloroethane', ['ethyl chloride'], 'CCCl', [0, 1], 'Halides'),
    P('1-chloropropane', '1-chloropropane', [], 'ClCCC', [1, 2], 'Halides'),
    P('1,2-dichloroethane', '1,2-dichloroethane', [], 'ClCCCl', [1, 2], 'Halides'),
    P('1,2-dibromoethane', '1,2-dibromoethane', [], 'BrCCBr', [1, 2], 'Halides'),
    P('1,2-difluoroethane', '1,2-difluoroethane', [], 'FCCF', [1, 2], 'Halides'),
    P('ethanol', 'ethanol', ['ethyl alcohol'], 'CCO', [0, 1], 'Alcohols & amines'),
    P('1-propanol', '1-propanol', ['propan-1-ol', 'n-propanol'], 'CCCO', [1, 2], 'Alcohols & amines'),
    P('2-butanol', '(R)-2-butanol', ['butan-2-ol', 'sec-butanol'], 'C[C@@H](O)CC', [1, 3], 'Alcohols & amines'),
    P('ethylene-glycol', 'ethylene glycol', ['1,2-ethanediol', 'ethane-1,2-diol'], 'OCCO', [1, 2], 'Alcohols & amines'),
    P('2-fluoroethanol', '2-fluoroethanol', [], 'FCCO', [1, 2], 'Alcohols & amines'),
    P('2-bromobutane', '(R)-2-bromobutane', [], 'C[C@@H](Br)CC', [1, 3], 'Stereo & E2'),
    P('meso-dibromobutane', 'meso-2,3-dibromobutane', ['(2R,3S)-2,3-dibromobutane'], 'C[C@@H](Br)[C@@H](Br)C', [1, 3], 'Stereo & E2'),
    P('rr-dibromobutane', '(2R,3R)-2,3-dibromobutane', [], 'C[C@@H](Br)[C@H](Br)C', [1, 3], 'Stereo & E2'),
    P('butylamine', '1-butanamine', ['butylamine'], 'NCCCC', [2, 3], 'Alcohols & amines'),
    P('propene', 'propene', ['propylene'], 'C=CC', [1, 2], 'Alkenes'),
    P('1-butene', '1-butene', ['but-1-ene'], 'C=CCC', [1, 2], 'Alkenes'),
    P('cyclohexane', 'cyclohexane', [], 'C1CCCCC1', [0, 1], 'Rings'),
    P('methylcyclohexane', 'methylcyclohexane', [], 'CC1CCCCC1', [1, 2], 'Rings'),
    P('cyclohexanol', 'cyclohexanol', [], 'OC1CCCCC1', [1, 2], 'Rings')
  ];

  /* Extra names accepted when typed (not shown as presets). [display name, SMILES, ...aliases] */
  const EXTRA_NAMES = [
    ['methane', 'C'], ['methanol', 'CO', 'methyl alcohol'],
    ['2,2-dimethylpropane', 'CC(C)(C)C', 'neopentane'],
    ['2-methylpentane', 'CC(C)CCC', 'isohexane'], ['3-methylpentane', 'CCC(C)CC'],
    ['heptane', 'CCCCCCC'], ['octane', 'CCCCCCCC'], ['nonane', 'CCCCCCCCC'], ['decane', 'CCCCCCCCCC'],
    ['2-propanol', 'CC(O)C', 'isopropanol', 'isopropyl alcohol', 'propan-2-ol'],
    ['1-butanol', 'CCCCO', 'butanol', 'butan-1-ol'],
    ['2-methyl-1-propanol', 'CC(C)CO', 'isobutanol', 'isobutyl alcohol'],
    ['2-methyl-2-propanol', 'CC(C)(C)O', 'tert-butanol', 't-butanol', 'tert-butyl alcohol'],
    ['1-pentanol', 'CCCCCO', 'pentanol'], ['1-hexanol', 'CCCCCCO', 'hexanol'],
    ['propanol', 'CCCO'],
    ['dimethyl ether', 'COC', 'methoxymethane'], ['diethyl ether', 'CCOCC', 'ether', 'ethoxyethane'],
    ['methoxyethane', 'CCOC', 'ethyl methyl ether'],
    ['1,1-dichloroethane', 'CC(Cl)Cl'], ['1,1-dibromoethane', 'CC(Br)Br'], ['1,1-difluoroethane', 'CC(F)F'],
    ['1-bromobutane', 'BrCCCC', 'butyl bromide'], ['1-bromo-2-methylpropane', 'BrCC(C)C', 'isobutyl bromide'],
    ['2-bromo-2-methylpropane', 'CC(C)(C)Br', 'tert-butyl bromide', 't-butyl bromide'],
    ['bromoethane', 'CCBr', 'ethyl bromide'], ['fluoroethane', 'CCF'], ['iodoethane', 'CCI', 'ethyl iodide'],
    ['chloromethane', 'CCl', 'methyl chloride'], ['2-chloropropane', 'CC(C)Cl', 'isopropyl chloride'],
    ['1-chlorobutane', 'ClCCCC'], ['2-chlorobutane', 'CC(Cl)CC'], ['1-bromopropane', 'BrCCC'],
    ['2-bromopropane', 'CC(C)Br'], ['1,2-dichloropropane', 'CC(Cl)CCl'],
    ['1,2-dibromobutane', 'BrCC(Br)CC'],
    ['2-chloroethanol', 'OCCCl'], ['2-aminoethanol', 'NCCO', 'ethanolamine'],
    ['ethylenediamine', 'NCCN', '1,2-diaminoethane'], ['2-methoxyethanol', 'COCCO'],
    ['1,2-dimethoxyethane', 'COCCOC', 'dme'], ['1,2-propanediol', 'CC(O)CO', 'propylene glycol'],
    ['methylamine', 'CN'], ['ethylamine', 'CCN', 'ethanamine'], ['propylamine', 'CCCN', '1-propanamine'],
    ['2-butanamine', 'CC(N)CC', 'sec-butylamine'], ['dimethylamine', 'CNC'], ['diethylamine', 'CCNCC'],
    ['trimethylamine', 'CN(C)C'], ['ethanethiol', 'CCS'], ['dimethyl sulfide', 'CSC'],
    ['ethene', 'C=C', 'ethylene'], ['ethyne', 'C#C', 'acetylene'], ['propyne', 'CC#C'],
    ['1-butyne', 'C#CCC'], ['2-butyne', 'CC#CC'], ['2-butene', 'CC=CC', 'but-2-ene'],
    ['trans-2-butene', 'C/C=C/C', '(E)-2-butene', 'E-2-butene', '(E)-but-2-ene'],
    ['cis-2-butene', 'C/C=C\\C', '(Z)-2-butene', 'Z-2-butene', '(Z)-but-2-ene'],
    ['trans-1,2-dichloroethene', 'Cl/C=C/Cl', '(E)-1,2-dichloroethene'], ['cis-1,2-dichloroethene', 'Cl/C=C\\Cl', '(Z)-1,2-dichloroethene'],
    ['2-methylpropene', 'C=C(C)C', 'isobutylene', 'isobutene'], ['1-pentene', 'C=CCCC'],
    ['1-hexene', 'C=CCCCC'], ['1-heptene', 'C=CCCCCC'], ['1,3-butadiene', 'C=CC=C', 'butadiene'],
    ['cyclopropane', 'C1CC1'], ['cyclobutane', 'C1CCC1'], ['cyclopentane', 'C1CCCC1'],
    ['methylcyclopentane', 'CC1CCCC1'], ['cycloheptane', 'C1CCCCCC1'], ['cyclohexene', 'C1=CCCCC1'],
    ['ethylcyclohexane', 'CCC1CCCCC1'], ['tert-butylcyclohexane', 'CC(C)(C)C1CCCCC1'],
    ['chlorocyclohexane', 'ClC1CCCCC1'], ['bromocyclohexane', 'BrC1CCCCC1'],
    ['cyclohexanone', 'O=C1CCCCC1'], ['cyclohexylamine', 'NC1CCCCC1'],
    ['benzene', 'c1ccccc1'], ['toluene', 'Cc1ccccc1', 'methylbenzene'], ['ethylbenzene', 'CCc1ccccc1'],
    ['propylbenzene', 'CCCc1ccccc1'], ['phenol', 'Oc1ccccc1'],
    ['acetaldehyde', 'CC=O', 'ethanal'], ['propanal', 'CCC=O', 'propionaldehyde'], ['butanal', 'CCCC=O'],
    ['acetone', 'CC(C)=O', 'propanone'], ['2-butanone', 'CCC(C)=O', 'butanone', 'methyl ethyl ketone'],
    ['acetic acid', 'CC(=O)O', 'ethanoic acid'], ['propanoic acid', 'CCC(=O)O', 'propionic acid'],
    ['butanoic acid', 'CCCC(=O)O', 'butyric acid'], ['methyl acetate', 'CC(=O)OC'],
    ['ethyl acetate', 'CCOC(C)=O'], ['acetonitrile', 'CC#N'], ['propanenitrile', 'CCC#N', 'propionitrile'],
    ['acetamide', 'CC(N)=O'], ['ethylene oxide', 'C1CO1', 'oxirane']
  ];
  /* Names that mean more than one molecule: return choices. */
  const AMBIGUOUS_NAMES = {
    '2,3dibromobutane': [
      { label: 'meso-2,3-dibromobutane (2R,3S)', value: 'meso-dibromobutane' },
      { label: '(2R,3R)-2,3-dibromobutane', value: 'rr-dibromobutane' }
    ]
  };

  function normName(s) {
    return String(s).toLowerCase().trim().replace(/\s+/g, ' ').replace(/^n[-\s]/, '').replace(/[\s\-_–]/g, '');
  }
  // default bond for typed names without a preset (C2=C3 in 2-butene: the fixed bond teaches)
  const EXTRA_DEFAULT = { 'CC=CC': [1, 2], 'C/C=C/C': [1, 2], 'C/C=C\\C': [1, 2] };
  const NAME_INDEX = new Map();   // normalized name -> {preset} | {name, smiles}
  const SMILES_INDEX = new Map(); // exact SMILES -> {name, presetId, defaultBond}
  (function buildIndexes() {
    const put = (k, v) => { const n = normName(k); if (n && !NAME_INDEX.has(n)) NAME_INDEX.set(n, v); };
    for (const p of presets) {
      const v = { preset: p };
      put(p.id, v); put(p.name, v);
      p.aliases.forEach(a => put(a, v));
      const bare = p.name.replace(/^\([RS]\)-/, '');
      if (bare !== p.name) put(bare, v);
      if (!SMILES_INDEX.has(p.smiles)) SMILES_INDEX.set(p.smiles, { name: p.name, presetId: p.id, defaultBond: p.defaultBond });
    }
    for (const row of EXTRA_NAMES) {
      const v = { name: row[0], smiles: row[1], defaultBond: EXTRA_DEFAULT[row[1]] || null };
      for (let k = 0; k < row.length; k++) if (k !== 1) put(row[k], v);
      if (!SMILES_INDEX.has(row[1])) SMILES_INDEX.set(row[1], { name: row[0], presetId: null, defaultBond: EXTRA_DEFAULT[row[1]] || null });
    }
  })();

  /* ------------------------------------------------------------------ */
  /* Molecular formula -> isomer choices                                  */
  /* ------------------------------------------------------------------ */
  function parseFormulaCounts(f) {
    const re = /([A-Z][a-z]?)(\d*)/g; const c = {}; let m;
    while ((m = re.exec(f))) { if (!m[1]) break; c[m[1]] = (c[m[1]] || 0) + (m[2] ? +m[2] : 1); }
    return c;
  }
  const ISOMERS = {
    CH4: [['methane', 'C']],
    C2H6: [['ethane', 'ethane']],
    C3H8: [['propane', 'propane']],
    C4H10: [['butane', 'butane'], ['2-methylpropane', 'isobutane']],
    C5H12: [['pentane', 'pentane'], ['2-methylbutane', '2-methylbutane'], ['2,2-dimethylpropane', 'CC(C)(C)C']],
    C6H14: [['hexane', 'hexane'], ['2-methylpentane', 'CC(C)CCC'], ['3-methylpentane', 'CCC(C)CC'],
      ['2,2-dimethylbutane', '2,2-dimethylbutane'], ['2,3-dimethylbutane', '2,3-dimethylbutane']],
    C2H6O: [['ethanol', 'ethanol'], ['dimethyl ether', 'COC']],
    C3H8O: [['1-propanol', '1-propanol'], ['2-propanol', 'CC(O)C'], ['methoxyethane', 'CCOC']],
    C4H10O: [['1-butanol', 'CCCCO'], ['2-butanol', '2-butanol'], ['2-methyl-1-propanol', 'CC(C)CO'],
      ['2-methyl-2-propanol', 'CC(C)(C)O'], ['diethyl ether', 'CCOCC']],
    C2H5Cl: [['chloroethane', 'chloroethane']],
    C2H5Br: [['bromoethane', 'CCBr']],
    C3H7Cl: [['1-chloropropane', '1-chloropropane'], ['2-chloropropane', 'CC(C)Cl']],
    C2H4Cl2: [['1,2-dichloroethane', '1,2-dichloroethane'], ['1,1-dichloroethane', 'CC(Cl)Cl']],
    C2H4Br2: [['1,2-dibromoethane', '1,2-dibromoethane'], ['1,1-dibromoethane', 'CC(Br)Br']],
    C2H4F2: [['1,2-difluoroethane', '1,2-difluoroethane'], ['1,1-difluoroethane', 'CC(F)F']],
    C4H9Br: [['1-bromobutane', 'BrCCCC'], ['2-bromobutane', '2-bromobutane'],
      ['1-bromo-2-methylpropane', 'BrCC(C)C'], ['2-bromo-2-methylpropane', 'CC(C)(C)Br']],
    C4H8Br2: [['meso-2,3-dibromobutane', 'meso-dibromobutane'], ['(2R,3R)-2,3-dibromobutane', 'rr-dibromobutane'],
      ['1,2-dibromobutane', 'BrCC(Br)CC']],
    C2H6O2: [['ethylene glycol', 'ethylene-glycol']],
    C2H5FO: [['2-fluoroethanol', '2-fluoroethanol']],
    C4H11N: [['1-butanamine', 'butylamine'], ['2-butanamine', 'CC(N)CC'], ['diethylamine', 'CCNCC']],
    C6H12: [['cyclohexane', 'cyclohexane'], ['1-hexene', 'C=CCCCC'], ['methylcyclopentane', 'CC1CCCC1']],
    C7H14: [['methylcyclohexane', 'methylcyclohexane'], ['1-heptene', 'C=CCCCCC']],
    C6H12O: [['cyclohexanol', 'cyclohexanol']],
    C3H6: [['propene', 'propene'], ['cyclopropane', 'C1CC1']],
    C4H8: [['1-butene', '1-butene'], ['2-butene', 'CC=CC'], ['2-methylpropene', 'C=C(C)C'], ['cyclobutane', 'C1CCC1']],
    C2H4: [['ethene', 'C=C']],
    C2H2: [['ethyne', 'C#C']],
    C6H6: [['benzene', 'c1ccccc1']],
    C7H8: [['toluene', 'Cc1ccccc1']]
  };
  const ALKANE_ROOT = { 5: 'pent', 6: 'hex', 7: 'hept', 8: 'oct', 9: 'non' };
  function isomersForFormula(f) {
    if (!f || !/^([A-Z][a-z]?\d*)+$/.test(f)) return null;
    const counts = parseFormulaCounts(f);
    const key = hill(counts);
    let rows = ISOMERS[key];
    if (!rows) {
      const els = Object.keys(counts);
      const n = counts.C || 0;
      // CnH2n+2 fallback: unbranched + single-methyl isomers
      if (els.length === 2 && n >= 7 && n <= 10 && counts.H === 2 * n + 2) {
        const names = { 7: 'heptane', 8: 'octane', 9: 'nonane', 10: 'decane' };
        rows = [[names[n], 'C'.repeat(n)]];
        const L = n - 1; // parent chain length for a methyl branch
        for (let pos = 2; pos <= Math.floor((L + 1) / 2) && rows.length < 6; pos++) {
          const smi = 'C'.repeat(pos - 1) + 'C(C)' + 'C'.repeat(L - pos);
          rows.push([pos + '-methyl' + ALKANE_ROOT[L] + 'ane', smi]);
        }
      }
    }
    if (!rows) return null;
    return rows.map(([name, value]) => {
      const p = presets.find(x => x.id === value);
      const smi = p ? p.smiles : value;
      const display = p ? p.name : name;
      return { label: display + ' · ' + smi, value };
    });
  }

  /* ------------------------------------------------------------------ */
  /* Shared finishing: H counts, validation, Mol assembly                 */
  /* ------------------------------------------------------------------ */
  function fail(error, extra) { return Object.assign({ ok: false, error }, extra || {}); }
  function elName(el) { return EL_NAME[el] || el; }
  function overMsg(el, idx, total, suffix) {
    const max = maxValence(el, 0);
    return elName(el) + ' ' + (idx + 1) + ' would have ' + total + ' bonds, but ' +
      elName(el).toLowerCase() + ' makes ' + max + '.' + (suffix || '');
  }

  /*
   * Fold explicit H atoms into their heavy neighbour. atoms may carry
   * `hCount` (explicit count) or null (to be computed), and `extraH` is added.
   * Returns {atoms, bonds} with heavy atoms only, or {error}.
   */
  function foldHydrogens(atoms, bonds) {
    const isH = atoms.map(a => a.el === 'H');
    if (!isH.some(Boolean)) return { atoms, bonds };
    const nb = atoms.map(() => []);
    bonds.forEach(b => { nb[b.a].push(b.b); nb[b.b].push(b.a); });
    for (let i = 0; i < atoms.length; i++) {
      if (!isH[i]) continue;
      const hs = nb[i];
      if (hs.length !== 1 || isH[hs[0]] || atoms[i].charge) {
        return { error: 'Every H has to be attached to exactly one heavy atom (C, N, O, ...).' };
      }
      const p = atoms[hs[0]];
      p.extraH = (p.extraH || 0) + 1;
      if (p.chiralNbrs) p.chiralNbrs = p.chiralNbrs.map(x => (x === i ? -1 : x));
    }
    const map = new Map(); const out = [];
    atoms.forEach((a, i) => { if (!isH[i]) { map.set(i, out.length); out.push(a); } });
    if (!out.length) return { error: 'Add at least one carbon or other heavy atom.' };
    out.forEach(a => {
      if (a.chiralNbrs) a.chiralNbrs = a.chiralNbrs.map(x => (x < 0 ? -1 : map.has(x) ? map.get(x) : -1));
    });
    const nbonds = bonds.filter(b => !isH[b.a] && !isH[b.b]).map(b => Object.assign({}, b, { a: map.get(b.a), b: map.get(b.b) }));
    return { atoms: out, bonds: nbonds };
  }

  /*
   * Compute hCount for atoms whose hCount is null (from valence), add extraH,
   * and check nobody is over valence. `suffix` is appended to error messages.
   */
  function assignHydrogens(atoms, bonds, suffix) {
    const sums = orderSums(atoms.length, bonds);
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      const extra = a.extraH || 0;
      const used = sums[i] + extra;
      if (a.hCount == null) {
        const v = fitValence(a.el, a.charge, used);
        if (v < 0) return fail(overMsg(a.el, i, used, suffix), { atom: i, kind: 'valence' });
        a.hCount = v - sums[i];
      } else {
        a.hCount += extra;
        const total = sums[i] + a.hCount;
        if (total > maxValence(a.el, a.charge)) {
          return fail(overMsg(a.el, i, total, suffix), { atom: i, kind: 'valence' });
        }
      }
      delete a.extraH;
    }
    return null;
  }

  function suggestBond(mol) {
    const n = mol.nHeavy;
    const deg = i => heavyNeighbors(mol, i).length + (mol.atoms[i].hTotal || 0);
    const cands = [];
    mol.bonds.forEach(b => {
      if (b.a >= n || b.b >= n) return;
      if (deg(b.a) < 2 || deg(b.b) < 2) return;
      cands.push(b);
    });
    if (!cands.length) return null;
    const rot = cands.filter(b => b.order === 1 && !b.inRing &&
      hybrid(mol, b.a) !== 'sp' && hybrid(mol, b.b) !== 'sp');
    const pair = b => [Math.min(b.a, b.b), Math.max(b.a, b.b)];
    if (rot.length) {
      let best = null, bestKey = null;
      for (const b of rot) {
        const s1 = subtree(mol, b.a, b.b).length, s2 = subtree(mol, b.b, b.a).length;
        const cc = (mol.atoms[b.a].el === 'C' && mol.atoms[b.b].el === 'C') ? 1 : 0;
        const [lo, hi] = pair(b);
        const key = [Math.min(s1, s2), cc, -lo, -hi];
        let better = !best;
        if (!better) for (let k = 0; k < key.length; k++) { if (key[k] !== bestKey[k]) { better = key[k] > bestKey[k]; break; } }
        if (better) { best = b; bestKey = key; }
      }
      return pair(best);
    }
    const ringB = cands.filter(b => b.order === 1 && b.inRing);
    const pick = (ringB.length ? ringB : cands).map(pair).sort((x, y) => x[0] - y[0] || x[1] - y[1])[0];
    return pick;
  }

  /* Assemble, perceive and sanity-check a heavy-atom graph. */
  function finalize(atoms, bonds, meta) {
    const source = meta.source;
    if (!atoms.length) return fail(source === 'sketch' ? 'Draw a structure first.' : 'Type a name, formula, or SMILES.');
    if (atoms.length > HEAVY_MAX) return fail("That's bigger than this tool handles (30 heavy atoms max).");
    const pieces = countPieces(atoms.length, bonds);
    if (pieces > 1) {
      return fail(source === 'sketch'
        ? 'Connect everything into one molecule (you have ' + pieces + ' separate pieces).'
        : 'Draw one connected molecule.');
    }
    const mol = {
      atoms: atoms.map(a => {
        const o = { el: a.el, charge: a.charge | 0, hCount: a.hCount, hTotal: a.hCount };
        if (a.chiral) { o.chiral = a.chiral; o.chiralNbrs = (a.chiralNbrs || []).slice(); }
        if (a.x != null) o.x = a.x;
        if (a.y != null) o.y = a.y;
        return o;
      }),
      bonds: bonds.map(b => (b.stereo ? { a: b.a, b: b.b, order: b.order, inRing: false, stereo: Object.assign({}, b.stereo) }
        : { a: b.a, b: b.b, order: b.order, inRing: false })),
      nHeavy: atoms.length,
      explicitH: false,
      rings: [],
      name: meta.name || null,
      source,
      input: meta.input || '',
      defaultBond: meta.defaultBond ? meta.defaultBond.slice() : null,
      presetId: meta.presetId || null
    };
    perceive(mol);
    if (!mol.name) { const nm = nameForGraph(mol); if (nm) mol.name = nm; }
    if (!mol.defaultBond) mol.defaultBond = suggestBond(mol);
    return { ok: true, mol };
  }

  /* Graph-based name lookup (integration fix): CH3CH2CH2CH3, a drawn zigzag or CC(C)CC get
   * the same name as the typed name. Weisfeiler-Lehman hash of the heavy-atom graph (element,
   * charge, H count, bond order). Stereo presets are indexed under their name without the
   * stereo descriptor, since the input carries no stereo. */
  let GRAPH_INDEX = null, buildingGraphIndex = false;
  function graphKey(mol) {
    const h = [];
    for (let i = 0; i < mol.atoms.length; i++) { const a = mol.atoms[i]; if (a.el !== 'H') h.push(i); }
    const pos = new Map(); h.forEach((i, k) => pos.set(i, k));
    const adj = h.map(() => []);
    for (const b of mol.bonds) {
      if (!pos.has(b.a) || !pos.has(b.b)) continue;
      adj[pos.get(b.a)].push([pos.get(b.b), b.order]); adj[pos.get(b.b)].push([pos.get(b.a), b.order]);
    }
    let lab = h.map((i) => { const a = mol.atoms[i]; return a.el + ':' + (a.charge | 0) + ':' + (a.hTotal != null ? a.hTotal : a.hCount); });
    for (let it = 0; it < 5; it++) {
      lab = lab.map((l, k) => l + '(' + adj[k].map(([n, o]) => o + lab[n]).sort().join(',') + ')');
      // compress
      const u = Array.from(new Set(lab)).sort(); const m = new Map(u.map((x, j) => [x, 'L' + j + '_' + x.length]));
      lab = lab.map((x) => x.length > 60 ? m.get(x) + x.slice(0, 40) : x);
    }
    return h.length + '|' + lab.slice().sort().join('/');
  }
  function nameForGraph(mol) {
    if (buildingGraphIndex) return null;
    if (!GRAPH_INDEX) {
      buildingGraphIndex = true; GRAPH_INDEX = new Map();
      try {
        const add = (smi, name) => {
          if (/@/.test(smi)) { smi = smi.replace(/@+/g, ''); name = name.replace(/^(\([^)]*\)-|meso-)/, ''); }
          const r = parseSmiles(smi); if (!r.ok) return;
          const f = finalize(r.atoms.map(a => ({ el: a.el, charge: a.charge, hCount: a.hCount })), r.bonds, { source: 'smiles' });
          if (!f.ok) return;
          const k = graphKey(f.mol); if (!GRAPH_INDEX.has(k)) GRAPH_INDEX.set(k, name);
        };
        presets.forEach((p) => add(p.smiles, p.name));
        EXTRA_NAMES.forEach((row) => add(row[1], row[0]));
      } catch (e) { /* naming is a convenience only */ }
      buildingGraphIndex = false;
    }
    try { return GRAPH_INDEX.get(graphKey(mol)) || null; } catch (e) { return null; }
  }

  /* ------------------------------------------------------------------ */
  /* fromGraph (sketch and any caller with a ready graph)                 */
  /* ------------------------------------------------------------------ */
  function fromGraph(atomsIn, bondsIn, source, extra) {
    extra = extra || {};
    source = source || 'sketch';
    if (!atomsIn || !atomsIn.length) return fail(source === 'sketch' ? 'Draw a structure first.' : 'Type a name, formula, or SMILES.');
    const atoms = [];
    for (let i = 0; i < atomsIn.length; i++) {
      const a = atomsIn[i] || {};
      if (!EL[a.el]) return fail("I don't know the element '" + a.el + "'.", { atom: i });
      const o = { el: a.el, charge: a.charge | 0, hCount: a.hCount == null ? null : a.hCount | 0 };
      if (a.x != null) o.x = a.x;
      if (a.y != null) o.y = a.y;
      if (a.chiral) { o.chiral = a.chiral; o.chiralNbrs = (a.chiralNbrs || []).slice(); }
      atoms.push(o);
    }
    const bonds = []; const seen = new Set();
    for (const b of bondsIn || []) {
      const x = b.a | 0, y = b.b | 0, order = b.order == null ? 1 : b.order | 0;
      if (x < 0 || y < 0 || x >= atoms.length || y >= atoms.length || x === y) return fail('A bond has to join two different atoms.');
      if (order < 1 || order > 3) return fail('Bonds can be single, double, or triple.');
      const key = Math.min(x, y) + '-' + Math.max(x, y);
      if (seen.has(key)) return fail('Two atoms are joined twice.');
      seen.add(key);
      bonds.push({ a: x, b: y, order });
    }
    const folded = foldHydrogens(atoms, bonds);
    if (folded.error) return fail(folded.error);
    const suffix = source === 'sketch' ? ' Remove a bond or change the atom.' : '';
    const err = assignHydrogens(folded.atoms, folded.bonds, suffix);
    if (err) return err;
    return finalize(folded.atoms, folded.bonds, {
      source, input: source === 'sketch' ? '' : (extra.input || ''),
      name: extra.name || null, defaultBond: extra.defaultBond || null, presetId: extra.presetId || null
    });
  }

  /* ------------------------------------------------------------------ */
  /* SMILES                                                               */
  /* ------------------------------------------------------------------ */
  const ORGANIC = { B: 1, C: 1, N: 1, O: 1, P: 1, S: 1, F: 1, I: 1, Cl: 1, Br: 1 };
  const AROMATIC = { b: 'B', c: 'C', n: 'N', o: 'O', p: 'P', s: 'S' };
  const SUPPORTED_TEXT = 'This tool handles H, B, C, N, O, P, S, F, Cl, Br, and I.';

  function parseSmiles(str) {
    const atoms = []; const bonds = [];
    const rings = {}; const stack = [];
    let prev = -1; let pending = null; let pendDir = null; let i = 0;
    const n = str.length;
    const near = k => "I can't read the SMILES near '" + str.slice(k, k + 4) + "' (character " + (k + 1) + ').';
    const bondExists = (x, y) => bonds.some(b => (b.a === x && b.b === y) || (b.a === y && b.b === x));
    const mkBond = (x, y, sym) => {
      let order = 1, arom = false;
      if (sym === '=') order = 2;
      else if (sym === '#') order = 3;
      else if (sym === ':') arom = true;
      else if (sym == null && atoms[x].aromatic && atoms[y].aromatic) arom = true;
      bonds.push({ a: x, b: y, order, arom, implicitArom: arom && sym == null });
    };
    while (i < n) {
      const ch = str[i];
      if (ch === '(') {
        if (prev < 0) return fail("A branch can't come before the first atom.");
        stack.push(prev); i++; continue;
      }
      if (ch === ')') {
        if (!stack.length) return fail("There's a ) without a matching (.");
        prev = stack.pop(); i++; continue;
      }
      if (ch === '-' || ch === '=' || ch === '#' || ch === ':' || ch === '/' || ch === '\\') {
        if (prev < 0) return fail(near(i));
        pending = (ch === '/' || ch === '\\') ? '-' : ch;
        pendDir = (ch === '/' || ch === '\\') ? ch : null;
        i++; continue;
      }
      if (ch === '$') return fail("Quadruple bonds aren't supported here.");
      if (ch === '.') return fail('Draw one connected molecule.');
      if ((ch >= '0' && ch <= '9') || ch === '%') {
        if (prev < 0) return fail(near(i));
        let num;
        if (ch === '%') {
          const d = str.slice(i + 1, i + 3);
          if (!/^\d\d$/.test(d)) return fail(near(i));
          num = +d; i += 3;
        } else { num = +ch; i++; }
        if (rings[num]) {
          const r = rings[num];
          delete rings[num];
          if (r.atom === prev) return fail('Ring number ' + num + ' connects an atom to itself.');
          if (bondExists(r.atom, prev)) return fail('Two atoms are joined twice (ring number ' + num + ').');
          const sym = pending || r.sym;
          if (pending && r.sym && pending !== r.sym) return fail('Ring number ' + num + ' has two different bond types.');
          mkBond(r.atom, prev, sym); pendDir = null;
          atoms[r.atom].nbr[r.slot] = prev;
          atoms[prev].nbr.push(r.atom);
        } else {
          rings[num] = { atom: prev, sym: pending, slot: atoms[prev].nbr.length };
          atoms[prev].nbr.push(null); // filled at closure
        }
        pending = null;
        continue;
      }
      // ---- an atom ----
      let atom = null;
      if (ch === '[') {
        const close = str.indexOf(']', i);
        if (close < 0) return fail('A [ bracket atom never closes with ].');
        const body = str.slice(i + 1, close);
        const m = /^(\d*)([A-Z][a-z]?|[a-z][a-z]?)(@{1,2}(?:TH[12]|AL[12]|SP[1-3]|TB\d{1,2}|OH\d{1,2})?)?(H\d*)?([+-]+\d*|[+-]\d+)?(:\d+)?$/.exec(body);
        if (!m) return fail("I couldn't read the bracket atom [" + body + '].');
        let sym = m[2]; let aromatic = false;
        if (/^[a-z]/.test(sym)) {
          if (!AROMATIC[sym]) return fail("I don't know the element '" + sym + "'. " + SUPPORTED_TEXT);
          sym = AROMATIC[sym]; aromatic = true;
        } else if (!EL[sym]) {
          // e.g. [CH3] would match 'C' + 'H3' only if regex split works; a two-letter unknown lands here
          return fail("I don't know the element '" + sym + "'. " + SUPPORTED_TEXT);
        }
        let charge = 0;
        if (m[5]) {
          const sgn = m[5][0] === '+' ? 1 : -1;
          const digits = m[5].replace(/[+-]/g, '');
          charge = sgn * (digits ? +digits : m[5].length);
        }
        const h = m[4] ? (m[4].length > 1 ? +m[4].slice(1) : 1) : 0;
        const chiral = m[3] ? (m[3].indexOf('@@') === 0 ? '@@' : '@') : null;
        atom = { el: sym, charge, hCount: h, aromatic, chiral, bracket: true, bracketH: h };
        i = close + 1;
      } else {
        const two = str.slice(i, i + 2);
        if (two === 'Cl' || two === 'Br') { atom = { el: two }; i += 2; }
        else if (ORGANIC[ch]) { atom = { el: ch }; i++; }
        else if (AROMATIC[ch]) { atom = { el: AROMATIC[ch], aromatic: true }; i++; }
        else if (/[A-Z]/.test(ch)) {
          const sym = /[a-z]/.test(str[i + 1] || '') ? ch + str[i + 1] : ch;
          return fail("I don't know the element '" + sym + "'.");
        } else if (i > 0 && /[A-Z]/.test(str[i - 1]) && /[a-z]/.test(ch)) {
          return fail("I don't know the element '" + str[i - 1] + ch + "'. " + SUPPORTED_TEXT);
        } else return fail(near(i));
        atom.charge = 0; atom.hCount = null; atom.chiral = null; atom.bracket = false;
      }
      const idx = atoms.length;
      atom.nbr = [];
      if (prev >= 0) atom.nbr.push(prev);
      if (atom.bracket && atom.bracketH === 1) atom.nbr.push(-1);
      atoms.push(atom);
      if (atoms.length > 4 * HEAVY_MAX) return fail("That's bigger than this tool handles (30 heavy atoms max).");
      if (prev >= 0) {
        mkBond(prev, idx, pending);
        if (pendDir) bonds[bonds.length - 1].dir = pendDir;   // written prev -> idx
        atoms[prev].nbr.push(idx);
      } else if (pending) return fail(near(i - 1));
      pending = null; pendDir = null;
      prev = idx;
    }
    if (stack.length) return fail('A branch opens with ( but never closes.');
    const open = Object.keys(rings);
    if (open.length) return fail('Unclosed ring: ring number ' + open[0] + ' opens but never closes.');
    if (pending) return fail('The SMILES ends with a bond symbol.');
    if (!atoms.length) return fail('Type a name, formula, or SMILES.');

    // chirality neighbour lists
    atoms.forEach(a => { if (a.chiral) a.chiralNbrs = a.nbr.map(x => (x == null ? -1 : x)); });

    // explicit [H] atoms -> hydrogen counts
    const folded = foldHydrogens(atoms, bonds);
    if (folded.error) return fail(folded.error);
    let A = folded.atoms, B = folded.bonds;
    if (A.length > HEAVY_MAX) return fail("That's bigger than this tool handles (30 heavy atoms max).");
    if (countPieces(A.length, B) > 1) return fail('Draw one connected molecule.');

    // aromatic bonds that are not in a ring are plain single bonds (e.g. biphenyl link)
    if (B.some(b => b.arom)) {
      const probe = A.map(a => ({ el: a.el }));
      const pb = B.map(b => ({ a: b.a, b: b.b, order: 1 }));
      perceiveRings(probe, pb);
      B.forEach((b, k) => { if (b.arom && b.implicitArom && !pb[k].inRing) b.arom = false; });
      const kerr = kekulize(A, B);
      if (kerr) return fail(kerr);
    }
    // cis/trans from / and \ marks: '/' written u -> v means v sits above u.
    B.forEach(db => {
      if (db.order !== 2) return;
      const side = (end, other) => {
        for (const q of B) {
          if (!q.dir || q === db) continue;
          let x = -1, rel = 0;
          if (q.a === end && q.b !== other) { x = q.b; rel = q.dir === '/' ? 1 : -1; }        // end written first: x above
          else if (q.b === end && q.a !== other) { x = q.a; rel = q.dir === '/' ? -1 : 1; }   // x written first: x below
          if (x >= 0) return { x, rel };
        }
        return null;
      };
      const sa = side(db.a, db.b), sb = side(db.b, db.a);
      if (sa && sb) db.stereo = { x: sa.x, y: sb.x, cis: sa.rel === sb.rel };
    });
    B = B.map(b => (b.stereo ? { a: b.a, b: b.b, order: b.order, stereo: b.stereo } : { a: b.a, b: b.b, order: b.order }));
    A = A.map(a => ({
      el: a.el, charge: a.charge | 0,
      hCount: a.bracket ? a.hCount : null, extraH: a.extraH || 0,
      chiral: a.chiral || null, chiralNbrs: a.chiral ? a.chiralNbrs : undefined
    }));
    const err = assignHydrogens(A, B, ' Check the bonds around it.');
    if (err) return err;
    return { ok: true, atoms: A, bonds: B };
  }

  /* Assign alternating double bonds to aromatic bonds (maximum matching with backtracking). */
  function kekulize(atoms, bonds) {
    const n = atoms.length;
    const sumNonArom = new Array(n).fill(0), aromDeg = new Array(n).fill(0), hasDouble = new Array(n).fill(false);
    bonds.forEach(b => {
      if (b.arom) { aromDeg[b.a]++; aromDeg[b.b]++; }
      else {
        sumNonArom[b.a] += b.order; sumNonArom[b.b] += b.order;
        if (b.order === 2) { hasDouble[b.a] = true; hasDouble[b.b] = true; }
      }
    });
    const need = new Array(n).fill(false);
    for (let i = 0; i < n; i++) {
      if (!aromDeg[i]) continue;
      const a = atoms[i];
      if (hasDouble[i]) continue;
      const used = sumNonArom[i] + aromDeg[i] + (a.bracket ? (a.hCount || 0) : 0) + (a.extraH || 0);
      const v = fitValence(a.el, a.charge, used);
      // an aromatic atom takes part in a double bond when it has a free valence left
      need[i] = v >= 0 ? (v - used) >= 1 : false;
    }
    const nbrs = atoms.map(() => []);
    bonds.forEach((b, k) => {
      if (b.arom && need[b.a] && need[b.b]) { nbrs[b.a].push({ n: b.b, k }); nbrs[b.b].push({ n: b.a, k }); }
    });
    const mate = new Array(n).fill(-1);
    const needList = [];
    for (let i = 0; i < n; i++) if (need[i]) needList.push(i);
    if (needList.length % 2) return "I couldn't draw alternating double bonds for that ring.";
    let steps = 0;
    function solve() {
      if (++steps > 200000) return false;
      let pick = -1, best = Infinity;
      for (const i of needList) {
        if (mate[i] >= 0) continue;
        let c = 0; for (const e of nbrs[i]) if (mate[e.n] < 0) c++;
        if (c < best) { best = c; pick = i; }
        if (c === 0) break;
      }
      if (pick < 0) return true;
      if (best === 0) return false;
      for (const e of nbrs[pick]) {
        if (mate[e.n] >= 0) continue;
        mate[pick] = e.k; mate[e.n] = e.k;
        if (solve()) return true;
        mate[pick] = -1; mate[e.n] = -1;
      }
      return false;
    }
    if (!solve()) return "I couldn't draw alternating double bonds for that ring.";
    bonds.forEach((b, k) => {
      if (!b.arom) return;
      b.order = (mate[b.a] === k && mate[b.b] === k) ? 2 : 1;
      b.arom = false;
    });
    return null;
  }

  /* ------------------------------------------------------------------ */
  /* Condensed formulas                                                   */
  /* ------------------------------------------------------------------ */
  function looksCondensed(s) {
    if (/[\[@%]/.test(s)) return false;
    return /(Cl|Br|[BCNOPSFI])H/.test(s) || /^H\d*(Cl|Br|[BCNOPSFI])/.test(s);
  }
  function rewriteCondensed(s) {
    return s.replace(/\s+/g, '')
      .replace(/tBu/g, 'C(CH3)3').replace(/iPr/g, 'CH(CH3)2')
      .replace(/Et/g, 'CH2CH3').replace(/Me/g, 'CH3')
      .replace(/C6H5/g, '§').replace(/Ph/g, '§')
      .replace(/C2H5/g, 'CH2CH3').replace(/CO2H/g, 'COOH').replace(/CO2/g, 'COO');
  }
  class CondErr extends Error {}

  function tokenizeCondensed(s) {
    const toks = []; let i = 0;
    const atStart = () => !toks.length || toks[toks.length - 1].t === '(' || toks[toks.length - 1].t === 'bond';
    const readDigits = () => { let d = ''; while (i < s.length && s[i] >= '0' && s[i] <= '9') d += s[i++]; return d; };
    const readEl = () => {
      const two = s.slice(i, i + 2);
      if (two === 'Cl' || two === 'Br') { i += 2; return two; }
      const ch = s[i];
      if (/[A-Z]/.test(ch || '')) {
        if (/[a-z]/.test(s[i + 1] || '')) throw new CondErr("I don't know the element '" + ch + s[i + 1] + "'.");
        if ('BCNOPSFI'.indexOf(ch) < 0) throw new CondErr("I don't know the element '" + ch + "'.");
        i++; return ch;
      }
      throw new CondErr("I can't read the formula near '" + s.slice(i, i + 3) + "'.");
    };
    while (i < s.length) {
      const ch = s[i];
      if (ch === '(') { toks.push({ t: '(' }); i++; continue; }
      if (ch === ')') { i++; const d = readDigits(); toks.push({ t: ')', n: d ? +d : 1 }); continue; }
      if (ch === '=') { toks.push({ t: 'bond', order: 2 }); i++; continue; }
      if (ch === '#' || ch === '≡') { toks.push({ t: 'bond', order: 3 }); i++; continue; }
      if (ch === '-' || ch === '–') { i++; continue; }
      if (ch === '§') { toks.push({ t: 'unit', el: 'C', h: 0, ph: true }); i++; continue; }
      if (ch === 'H') {
        if (!atStart()) throw new CondErr("I can't read the H at character " + (i + 1) + '. Write H counts right after their atom, like CH3.');
        i++; const d = readDigits();
        const el = readEl();
        toks.push({ t: 'unit', el, h: d ? +d : 1 });
        continue;
      }
      const el = readEl();
      let h = 0;
      if (s[i] === 'H') { i++; const d = readDigits(); h = d ? +d : 1; }
      const cnt = readDigits();
      const reps = cnt ? +cnt : 1;
      if (reps > HEAVY_MAX) throw new CondErr("That's bigger than this tool handles (30 heavy atoms max).");
      for (let r = 0; r < reps; r++) toks.push({ t: 'unit', el, h });
    }
    // match parentheses
    const st = [];
    toks.forEach((t, k) => {
      if (t.t === '(') st.push(k);
      else if (t.t === ')') {
        if (!st.length) throw new CondErr("There's a ) without a matching (.");
        const o = st.pop(); toks[o].match = k;
      }
    });
    if (st.length) throw new CondErr('A ( never closes with ).');
    return expandInline(toks);
  }

  /* (CH2)n between backbone units -> n chain units. Only simple groups of divalent units. */
  function expandInline(toks) {
    const out = [];
    for (let k = 0; k < toks.length; k++) {
      const t = toks[k];
      if (t.t === '(') {
        const inner = toks.slice(k + 1, t.match);
        const simple = inner.length && inner.every(u => u.t === 'unit' && !u.ph &&
          valenceList(u.el, 0)[0] - u.h === 2);
        if (simple) {
          const reps = toks[t.match].n;
          for (let r = 0; r < reps; r++) inner.forEach(u => out.push(Object.assign({}, u)));
          k = t.match;
          continue;
        }
      }
      out.push(Object.assign({}, t));
    }
    // recompute matches
    const st = [];
    out.forEach((t, k) => { if (t.t === '(') st.push(k); else if (t.t === ')') out[st.pop()].match = k; });
    return out;
  }

  function buildCondensed(toks, mask) {
    const atoms = [], bonds = [];
    const addAtom = (el, h) => { atoms.push({ el, charge: 0, hCount: h }); return atoms.length - 1; };
    const sumOf = i => { let s = 0; for (const b of bonds) if (b.a === i || b.b === i) s += b.order; return s; };
    const free = i => valenceList(atoms[i].el, 0)[0] - atoms[i].hCount - sumOf(i);
    const addBond = (a, b, o) => bonds.push({ a, b, order: o });
    const placeUnit = tk => {
      if (tk.ph) { // phenyl, Kekulé
        const r0 = addAtom('C', 0); const ring = [r0];
        for (let k = 0; k < 5; k++) ring.push(addAtom('C', 1));
        for (let k = 0; k < 6; k++) addBond(ring[k], ring[(k + 1) % 6], k % 2 === 0 ? 2 : 1);
        return r0;
      }
      return addAtom(tk.el, tk.h);
    };
    function seq(from, to, anchor) {
      const chain = []; let order = 1; let last = null; let deferred = [];
      let p = from;
      while (p < to) {
        const tk = toks[p];
        if (tk.t === 'bond') { order = tk.order; p++; continue; }
        if (tk.t === '(') {
          const close = tk.match, reps = toks[close].n;
          if (last == null) deferred.push({ s: p + 1, e: close, reps });
          else for (let r = 0; r < reps; r++) seq(p + 1, close, last);
          p = close + 1; continue;
        }
        if (tk.t === ')') throw new CondErr("There's a ) without a matching (.");
        const idx = placeUnit(tk);
        const terminal = tk.amb != null && ((mask >> tk.amb) & 1) === 1;
        if (last == null) {
          if (anchor != null) addBond(anchor, idx, order);
        } else {
          let target = -1;
          for (let k = chain.length - 1; k >= 0; k--) if (free(chain[k]) > 0) { target = chain[k]; break; }
          if (target < 0) target = chain.length ? chain[chain.length - 1] : last;
          addBond(target, idx, terminal ? 2 : order);
        }
        order = 1; last = idx;
        if (!terminal) chain.push(idx);
        for (const d of deferred) for (let r = 0; r < d.reps; r++) seq(d.s, d.e, idx);
        deferred = [];
        p++;
      }
      if (deferred.length) throw new CondErr('A group in parentheses has nothing to attach to.');
    }
    seq(0, toks.length, null);
    if (atoms.length > HEAVY_MAX) throw new CondErr("That's bigger than this tool handles (30 heavy atoms max).");

    // raise bond orders between atoms that are both short of valence (CH2CHCH3 -> propene)
    const sorted = bonds.slice().sort((x, y) => Math.min(x.a, x.b) - Math.min(y.a, y.b) || Math.max(x.a, x.b) - Math.max(y.a, y.b));
    for (const b of sorted) {
      while (b.order < 3 && free(b.a) > 0 && free(b.b) > 0) b.order++;
    }
    // validate: every atom must land exactly on an allowed valence
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      const total = sumOf(i) + a.hCount;
      const list = valenceList(a.el, 0);
      if (total > list[list.length - 1]) {
        return fail(overMsg(a.el, i, total, ' Check the H counts.'), { atom: i, kind: 'valence' });
      }
      if (list.indexOf(total) < 0) {
        const want = list.find(v => v > total);
        return fail(elName(a.el) + ' ' + (i + 1) + ' needs ' + want + ' bonds but has ' + total + '. Check the H counts.', { atom: i, kind: 'valence' });
      }
    }
    return { ok: true, atoms, bonds };
  }

  function parseCondensed(input) {
    let toks;
    try { toks = tokenizeCondensed(rewriteCondensed(input)); } catch (e) {
      if (e instanceof CondErr) return fail(e.message, { kind: 'syntax' });
      throw e;
    }
    if (!toks.some(t => t.t === 'unit')) return fail("I can't read that formula.", { kind: 'syntax' });
    // bare O/S directly after a unit (and not the final unit) could be a chain O or a C=O
    let lastUnit = -1; toks.forEach((t, k) => { if (t.t === 'unit') lastUnit = k; });
    let nAmb = 0;
    toks.forEach((t, k) => {
      if (t.t === 'unit' && !t.ph && (t.el === 'O' || t.el === 'S') && t.h === 0 && k !== lastUnit &&
        k > 0 && toks[k - 1].t === 'unit' && nAmb < 10) t.amb = nAmb++;
    });
    const masks = [];
    for (let m = 0; m < (1 << nAmb); m++) masks.push(m);
    const pop = m => { let c = 0; while (m) { c += m & 1; m >>= 1; } return c; };
    masks.sort((x, y) => pop(x) - pop(y) || x - y);
    let firstErr = null;
    for (const m of masks) {
      let r;
      try { r = buildCondensed(toks, m); } catch (e) {
        if (e instanceof CondErr) return fail(e.message, { kind: 'syntax' });
        throw e;
      }
      if (r.ok) return r;
      if (!firstErr) firstErr = r;
    }
    return firstErr;
  }

  /* ------------------------------------------------------------------ */
  /* parse                                                                */
  /* ------------------------------------------------------------------ */
  function fromParts(r, meta) {
    const A = r.atoms.map(a => ({ el: a.el, charge: a.charge, hCount: a.hCount, chiral: a.chiral, chiralNbrs: a.chiralNbrs }));
    return finalize(A, r.bonds, meta);
  }
  function parseSmilesMol(smiles, meta) {
    const r = parseSmiles(smiles);
    if (!r.ok) return r;
    return fromParts(r, meta);
  }

  function parse(input) {
    const raw = input == null ? '' : String(input);
    const s = raw.trim().replace(/\s+/g, ' ');
    if (!s) return fail('Type a name, formula, or SMILES.');

    // 1. preset / known name
    const key = normName(s);
    if (AMBIGUOUS_NAMES[key]) {
      return { ok: false, choices: AMBIGUOUS_NAMES[key].map(c => Object.assign({}, c)), message: s + ' could be more than one molecule. Which one?' };
    }
    const hit = NAME_INDEX.get(key);
    if (hit) {
      if (hit.preset) {
        const p = hit.preset;
        return parseSmilesMol(p.smiles, { source: 'name', input: raw, name: p.name, presetId: p.id, defaultBond: p.defaultBond });
      }
      return parseSmilesMol(hit.smiles, { source: 'name', input: raw, name: hit.name, defaultBond: hit.defaultBond || null });
    }

    const compact = s.replace(/\s+/g, '');
    // 2. molecular formula (each element once, contains H)
    if (/^([A-Z][a-z]?\d*)+$/.test(compact)) {
      const syms = compact.match(/[A-Z][a-z]?/g);
      const unique = new Set(syms).size === syms.length;
      if (unique && syms.indexOf('H') >= 0) {
        const bad = syms.find(x => !EL[x]);
        const choices = bad ? null : isomersForFormula(compact);
        if (bad) { /* e.g. EtOH: not a formula, let the condensed parser try */ }
        else if (choices && choices.length === 1) {
          const r = parse(choices[0].value);
          if (r.ok) { r.mol.source = 'formula'; r.mol.input = raw; }
          return r;
        }
        else if (choices && choices.length > 1) {
          return { ok: false, choices, message: compact + ' fits more than one molecule. Which one?' };
        }
        else {
          const c = parseCondensed(compact);
          if (c.ok) return fromParts(c, { source: 'condensed', input: raw });
          if (c.kind === 'valence') return fail(c.error, { atom: c.atom });
          return fail("I can't tell which molecule " + compact + " is. Type its name, a condensed formula like CH3CH2OH, or draw it.");
        }
      }
    }

    // 3. condensed formula
    if (looksCondensed(compact)) {
      const c = parseCondensed(compact);
      if (c.ok) return fromParts(c, { source: 'condensed', input: raw });
      return fail(c.error, c.atom != null ? { atom: c.atom } : null);
    }

    // 4. SMILES
    const known = SMILES_INDEX.get(compact);
    const r = parseSmiles(compact);
    if (r.ok) {
      return fromParts(r, {
        source: 'smiles', input: raw,
        name: known ? known.name : null,
        presetId: known ? known.presetId : null,
        defaultBond: known ? known.defaultBond : null
      });
    }
    // SMILES failed: maybe a condensed formula without H (EtCl, CCl3CH3) or an unknown name
    if (/^[A-Za-z0-9()=#≡\-]+$/.test(compact)) {
      const c = parseCondensed(compact);
      if (c.ok) return fromParts(c, { source: 'condensed', input: raw });
    }
    if (/[a-z]{3,}/.test(s) && /[adefghjkmqtuvwxyz]/.test(s.replace(/Cl|Br/g, ''))) {
      return fail("I don't know the name “" + s + '” yet. Try a condensed formula (CH3CH2CH2OH), a SMILES string (CCCO), or pick a preset.');
    }
    return fail(r.error, r.atom != null ? { atom: r.atom } : null);
  }

  /* ------------------------------------------------------------------ */
  const api = {
    EL,
    presets,
    parse,
    fromGraph,
    perceive,
    addHydrogens,
    neighbors,
    bondBetween,
    formula,
    atomLabel,
    bondLabel,
    groupInfo,
    groupLabel,
    hybrid,
    isomersForFormula,
    // extras (not in the spec contract, safe to ignore)
    valenceList,
    elementName: elName
  };
  root.NNChem = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

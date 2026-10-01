/* rs/cip.js: window.RSCip, the CIP (R/S) engine for the R/S Assigner.
 * Pure, no DOM, node-requirable. Contract: rs/SPEC.md section 3.
 * Every label it produces was checked against RDKit's new CIP labeler (rs/test/cip.test.js).
 *
 * Deviations from / additions to SPEC.md (documented here per the Newman convention):
 *  - classifyAtom checks 'tooComplex' (a capped comparison) before 'pseudo'/'sameGroups'/'ringSame', because a capped
 *    comparison returns 0 and would otherwise be misreported as a tie.
 *  - Classification messages come back with {list} and {A}/{B} already filled (group labels from NNChem.groupLabel,
 *    "ring " prefix when the group and the center share a ring). Only {C} (the center name) is left for the page.
 *    A page that still runs .replace('{A}', ...) is harmless (no-op).
 *  - 'pseudo' also covers ties whose two branches carry stereocenters that are neither identical nor mirror images
 *    (diastereomorphic branches need CIP Rule 4, out of scope).
 *  - Extra fields (not in the contract, safe to ignore): SetBox.fromEl, SetBox.terminal (true for boxes under H,
 *    duplicate or phantom nodes, which never hold atoms); Comparison.boxes = { a: SetBox, b: SetBox } (the deciding
 *    boxes); Ranking.ligands[i].ring; Classification.groups (labels in neighbor order).
 *  - Comparison.dupUsed = the WINNING deciding box contains a duplicate atom. The literal wording in 3.1 ("one of the two
 *    deciding boxes") contradicts the 3.7 table (cysteine: false although COOH's (O,O,O) holds a duplicate); this reading
 *    matches all 11 rows. Extra Comparison.dupDecisive = a duplicate at or before the first differing position.
 *  - For a tie, trailing rows that contain only boxes under H/duplicate atoms are dropped from Comparison.rows.
 *  - Extra function: explain(comparison, names?) -> string[], one plain sentence per sphere, for step-by-step display.
 *  - describe(): with several centers and no locants, text = base + ' (' + descriptor + ')' instead of descriptor-base.
 */
(function (root) {
  'use strict';

  const C = () => root.NNChem || require('../newman/chem.js');

  const Z = { H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, P: 15, S: 16, Cl: 17, Br: 35, I: 53 };
  const EN = { H: 2.20, B: 2.04, C: 2.55, N: 3.04, O: 3.44, F: 3.98, P: 2.19, S: 2.58, Cl: 3.16, Br: 2.96, I: 2.66 };
  const ZX = { Li: 3, Be: 4, Na: 11, Mg: 12, Al: 13, Si: 14, K: 19, Se: 34, Sn: 50 };
  const zOf = el => Z[el] || ZX[el] || 0;
  const ELNAME = { O: 'oxygen', N: 'nitrogen', F: 'fluorine', Cl: 'chlorine', Br: 'bromine', I: 'iodine', S: 'sulfur',
    P: 'phosphorus', B: 'boron', H: 'hydrogen' };

  const MAX_SPHERE = 12;
  const MAX_NODES = 5000;

  /* ------------------------------------------------------------------ */
  /* Graph helpers                                                        */
  /* ------------------------------------------------------------------ */
  const hCache = typeof WeakMap !== 'undefined' ? new WeakMap() : null;
  const adjCache = typeof WeakMap !== 'undefined' ? new WeakMap() : null;

  function withH(mol) {
    if (mol.explicitH) return mol;
    if (hCache && hCache.has(mol)) return hCache.get(mol);
    const m = C().addHydrogens(mol);
    if (hCache) hCache.set(mol, m);
    return m;
  }
  function adjOf(mol) {
    if (adjCache && adjCache.has(mol)) return adjCache.get(mol);
    const adj = mol.atoms.map(() => []);
    mol.bonds.forEach(b => {
      adj[b.a].push({ to: b.b, order: b.order || 1 });
      adj[b.b].push({ to: b.a, order: b.order || 1 });
    });
    adj.forEach(l => l.sort((x, y) => x.to - y.to));
    if (adjCache) adjCache.set(mol, adj);
    return adj;
  }
  const nbrs = (mol, i) => adjOf(mol)[i].map(e => e.to);
  const nHeavyOf = mol => (mol.nHeavy != null ? mol.nHeavy : mol.atoms.filter(a => a.el !== 'H').length);

  // atoms reachable from `start` without passing through `block`
  function subtree(mol, block, start) {
    const adj = adjOf(mol), seen = new Set([block, start]), out = [start], st = [start];
    while (st.length) {
      const x = st.pop();
      for (const e of adj[x]) if (!seen.has(e.to)) { seen.add(e.to); out.push(e.to); st.push(e.to); }
    }
    return out;
  }
  function shareRing(mol, a, b, c) {
    const rings = mol.rings || [];
    return rings.some(r => r.indexOf(a) >= 0 && r.indexOf(b) >= 0 && (c == null || r.indexOf(c) >= 0));
  }
  function groupName(mol, lig, center) {
    if (mol.atoms[lig].el === 'H') return 'H';
    let lab;
    try { lab = C().groupLabel(mol, lig, center); } catch (e) { lab = mol.atoms[lig].el; }
    if (shareRing(mol, lig, center)) lab = 'ring ' + lab;
    return lab;
  }
  function joinList(xs) {
    if (xs.length <= 1) return xs.join('');
    return xs.slice(0, -1).join(', ') + ' and ' + xs[xs.length - 1];
  }

  /* ------------------------------------------------------------------ */
  /* Hierarchical digraph (lazy) and the breadth-first comparison (3.1)   */
  /* ------------------------------------------------------------------ */
  function makeCtx(mol, opts) {
    return { mol, adj: adjOf(mol), dups: !(opts && opts.duplicates === false), nodes: 0, capped: false };
  }
  function newNode(ctx, atom, dup, path, parentAtom) {
    ctx.nodes++;
    if (ctx.nodes > MAX_NODES) ctx.capped = true;
    return { atom, el: ctx.mol.atoms[atom].el, z: zOf(ctx.mol.atoms[atom].el), dup, path, parentAtom, kids: null, sorted: null };
  }
  function rootNode(ctx, center, lig) { return newNode(ctx, lig, false, new Set([center, lig]), center); }

  function kidsOf(ctx, node) {
    if (node.kids) return node.kids;
    const out = [];
    node.kids = out;
    if (node.dup || node.z <= 1) return out;
    for (const e of ctx.adj[node.atom]) {
      const t = e.to;
      if (t === node.parentAtom) {
        if (ctx.dups) for (let k = 1; k < e.order; k++) out.push(newNode(ctx, t, true, null, node.atom));
        continue;
      }
      if (node.path.has(t)) {                       // ring closure
        if (ctx.dups) out.push(newNode(ctx, t, true, null, node.atom));
        continue;
      }
      const p = new Set(node.path); p.add(t);
      out.push(newNode(ctx, t, false, p, node.atom));
      if (ctx.dups) for (let k = 1; k < e.order; k++) out.push(newNode(ctx, t, true, null, node.atom));
    }
    return out;
  }
  function sortedKids(ctx, node) {
    if (node.sorted) return node.sorted;
    const k = kidsOf(ctx, node).slice();
    if (k.length > 1) k.sort((a, b) => cmpNodes(ctx, b, a, null));
    node.sorted = k;
    return k;
  }

  const entryOf = n => (n ? { el: n.el, atom: n.atom, dup: !!n.dup } : { el: '0', atom: null, dup: false });
  function boxOf(node, kids, len) {
    const atoms = [];
    for (let j = 0; j < len; j++) atoms.push(entryOf(kids[j]));
    return { from: node ? node.atom : null, fromEl: node ? node.el : '0',
      terminal: !node || node.dup || node.z <= 1, atoms };
  }

  // > 0 if A outranks B. With `tr` (an object) records rows / deciding info.
  function cmpNodes(ctx, A, B, tr) {
    if (ctx.capped) return 0;
    if (A.z !== B.z) {
      if (tr) {
        const ra = { from: null, fromEl: null, terminal: false, atoms: [entryOf(A)] };
        const rb = { from: null, fromEl: null, terminal: false, atoms: [entryOf(B)] };
        tr.rows.push({ sphere: 1, a: [ra], b: [rb], diff: { box: 0, pos: 0 } });
        tr.sphere = 1; tr.deciding = { a: ra.atoms[0], b: rb.atoms[0] }; tr.boxes = { a: ra, b: rb };
      }
      return A.z - B.z;
    }
    if (tr) {
      tr.rows.push({ sphere: 1, a: [{ from: null, fromEl: null, terminal: false, atoms: [entryOf(A)] }],
        b: [{ from: null, fromEl: null, terminal: false, atoms: [entryOf(B)] }], diff: null });
    }
    let LA = [A], LB = [B], sphere = 2;
    while (LA.length || LB.length) {
      if (sphere > MAX_SPHERE) { ctx.capped = true; return 0; }
      const n = Math.max(LA.length, LB.length);
      const row = tr ? { sphere, a: [], b: [], diff: null } : null;
      let res = 0;
      for (let i = 0; i < n; i++) {
        const ka = LA[i] ? sortedKids(ctx, LA[i]) : [];
        const kb = LB[i] ? sortedKids(ctx, LB[i]) : [];
        if (ctx.capped) return 0;
        const len = Math.max(3, ka.length, kb.length);
        if (row) { row.a.push(boxOf(LA[i], ka, len)); row.b.push(boxOf(LB[i], kb, len)); }
        if (res === 0) {
          for (let j = 0; j < len; j++) {
            const za = ka[j] ? ka[j].z : 0, zb = kb[j] ? kb[j].z : 0;
            if (za !== zb) {
              res = za - zb;
              if (row) {
                row.diff = { box: i, pos: j };
                tr.deciding = { a: row.a[i].atoms[j], b: row.b[i].atoms[j] };
                tr.boxes = { a: row.a[i], b: row.b[i] };
              }
              break;
            }
          }
          if (res !== 0 && !row) return res;
        }
      }
      if (row) tr.rows.push(row);
      if (res !== 0) { if (tr) tr.sphere = sphere; return res; }
      const NA = [], NB = [];
      for (const x of LA) NA.push.apply(NA, sortedKids(ctx, x));
      for (const x of LB) NB.push.apply(NB, sortedKids(ctx, x));
      LA = NA; LB = NB;
      sphere++;
    }
    return 0;
  }

  function setText(box) {
    if (!box || !box.atoms) return '()';
    return '(' + box.atoms.map(e => e.el).join(',') + ')';
  }

  function traced(ctx, A, B, a, b) {
    const tr = { rows: [], sphere: null, deciding: null, boxes: null };
    const s = cmpNodes(ctx, A, B, tr);
    const sign = ctx.capped ? 0 : Math.sign(s);
    let text = '';
    if (sign !== 0 && tr.boxes) {
      text = tr.sphere === 1 ? tr.deciding.a.el + ' vs ' + tr.deciding.b.el
        : setText(tr.boxes.a) + ' vs ' + setText(tr.boxes.b);
    } else {
      // drop trailing rows that hold no atoms at all (only boxes under H / duplicates)
      while (tr.rows.length > 1) {
        const r = tr.rows[tr.rows.length - 1];
        if (r.a.concat(r.b).every(bx => bx.terminal)) tr.rows.pop(); else break;
      }
      const last = tr.rows[tr.rows.length - 1];
      text = last ? (last.sphere === 1 ? A.el + ' vs ' + B.el
        : last.a.filter(bx => !bx.terminal).map(setText).join(' ') + ' vs ' +
          last.b.filter(bx => !bx.terminal).map(setText).join(' ')) : '';
    }
    // dupUsed: the winning deciding box holds a duplicate atom (matches every row of SPEC 3.7).
    // dupDecisive (extra): a duplicate sits at or before the first differing position in either deciding box.
    let dupUsed = false, dupDecisive = false;
    if (sign !== 0 && tr.boxes && tr.sphere > 1) {
      dupUsed = (sign > 0 ? tr.boxes.a : tr.boxes.b).atoms.some(e => e.dup);
      const last = tr.rows[tr.rows.length - 1], pos = last.diff ? last.diff.pos : 0;
      dupDecisive = tr.boxes.a.atoms.slice(0, pos + 1).some(e => e.dup) || tr.boxes.b.atoms.slice(0, pos + 1).some(e => e.dup);
    }
    return {
      a, b, sign,
      sphere: sign !== 0 ? tr.sphere : null,
      rows: tr.rows,
      text,
      deciding: sign !== 0 ? tr.deciding : null,
      boxes: sign !== 0 ? tr.boxes : null,
      dupUsed,
      dupDecisive,
      capped: ctx.capped
    };
  }

  function compareLigands(mol, center, a, b, opts) {
    const m = withH(mol);
    const ctx = makeCtx(m, opts);
    return traced(ctx, rootNode(ctx, center, a), rootNode(ctx, center, b), a, b);
  }

  /* ------------------------------------------------------------------ */
  /* Ranking                                                              */
  /* ------------------------------------------------------------------ */
  function rankSubstituents(mol, center, opts) {
    const m = withH(mol);
    const ctx = makeCtx(m, opts);
    const ligs = nbrs(m, center);
    const nodes = ligs.map(l => rootNode(ctx, center, l));
    const wins = ligs.map(() => 0);
    const ties = [];
    for (let i = 0; i < ligs.length; i++) {
      for (let j = i + 1; j < ligs.length; j++) {
        const s = cmpNodes(ctx, nodes[i], nodes[j], null);
        if (ctx.capped) break;
        if (s > 0) wins[i]++; else if (s < 0) wins[j]++; else ties.push([ligs[i], ligs[j]]);
      }
    }
    const capped = ctx.capped;
    const idx = ligs.map((_, i) => i).sort((x, y) => (wins[y] - wins[x]) || (x - y));
    const order = idx.map(i => ligs[i]);
    const tied = new Set();
    ties.forEach(p => { tied.add(p[0]); tied.add(p[1]); });
    const ligands = ligs.map(l => ({
      atom: l,
      rank: capped || tied.has(l) ? null : order.indexOf(l) + 1,
      el: m.atoms[l].el,
      label: groupName(m, l, center).replace(/^ring /, ''),
      ring: shareRing(m, l, center)
    }));
    const steps = [];
    for (let k = 0; k + 1 < order.length; k++) {
      const c2 = makeCtx(m, opts);
      steps.push(traced(c2, rootNode(c2, center, order[k]), rootNode(c2, center, order[k + 1]), order[k], order[k + 1]));
    }
    return { center, order, ligands, tie: capped || ties.length > 0, ties, capped, steps };
  }

  /* ------------------------------------------------------------------ */
  /* Chirality (3.3)                                                      */
  /* ------------------------------------------------------------------ */
  const sub = (p, q) => [p[0] - q[0], p[1] - q[1], (p[2] || 0) - (q[2] || 0)];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  function det3(a, b, c) {
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
  }
  function unit(v) { const l = Math.hypot(v[0], v[1], v[2]); return l > 0 ? [v[0] / l, v[1] / l, v[2] / l] : [0, 0, 0]; }

  function labelFromCoords(coords, order) {
    if (!coords || order.length !== 4) return null;
    const p = order.map(a => coords[a]);
    if (p.some(x => !x)) return null;
    const d = det3(sub(p[0], p[3]), sub(p[1], p[3]), sub(p[2], p[3]));
    if (Math.abs(d) < 1e-6) return null;
    return d < 0 ? 'R' : 'S';
  }

  function tagsOrder(m, center) {
    const a = m.atoms[center];
    if (!a || !a.chiral || !a.chiralNbrs || a.chiralNbrs.length !== 4) return null;
    const nb = a.chiralNbrs;
    if (nb.some(x => x == null || x < 0)) return null;
    const real = nbrs(m, center);
    if (real.length !== 4) return null;
    if (nb.slice().sort((x, y) => x - y).join() !== real.join()) return null;
    return nb;
  }
  function stereoGiven(mol, center) { return !!tagsOrder(withH(mol), center); }

  function tagLabel(m, center, order) {
    const nb = tagsOrder(m, center);
    if (!nb || order.length !== 4) return null;
    const perm = order.map(x => nb.indexOf(x));
    if (perm.some(p => p < 0)) return null;
    let inv = 0;
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) if (perm[i] > perm[j]) inv++;
    let at = m.atoms[center].chiral === '@';
    if (inv % 2) at = !at;          // permuted to rank order
    at = !at;                       // move r4 to the front: 3 transpositions
    return at ? 'R' : 'S';          // from r4, 1-2-3 CCW ('@') => with 4 away, clockwise => R
  }

  function chiralityFromTags(mol, center) {
    const m = withH(mol);
    const r = rankSubstituents(m, center);
    if (r.tie || r.order.length !== 4) return null;
    return tagLabel(m, center, r.order);
  }
  function chirality(mol, center, coords) {
    const m = withH(mol);
    const r = rankSubstituents(m, center);
    if (r.tie || r.order.length !== 4) return null;
    if (coords && coords.length >= m.atoms.length) return labelFromCoords(coords, r.order);
    return tagLabel(m, center, r.order);
  }

  /* ------------------------------------------------------------------ */
  /* Classification (3.2)                                                 */
  /* ------------------------------------------------------------------ */
  function cls(code, msg, extra) { return Object.assign({ code, msg }, extra || {}); }

  function classifyAtom(mol, i, coords, _busy) {
    const m = withH(mol);
    const a = m.atoms[i];
    if (!a) return cls('notCarbon', 'Only carbons count as stereocenters here.');
    const adj = adjOf(m)[i];
    if (a.el !== 'C') {
      if (a.el === 'O') return cls('notCarbon', "That's an oxygen. It has only two bonds, so it can't hold four different groups.");
      if (a.el === 'N' && !a.charge && adj.length === 3) {
        return cls('notCarbon', "That's a nitrogen. Its lone pair lets it flip inside out (invert) quickly, so we don't give it R or S.");
      }
      if (['F', 'Cl', 'Br', 'I'].indexOf(a.el) >= 0) return cls('notCarbon', "That's a " + ELNAME[a.el] + '. It has only one bond.');
      return cls('notCarbon', 'Only carbons count as stereocenters here.');
    }
    const maxOrder = adj.reduce((x, e) => Math.max(x, e.order), 1);
    if (adj.length < 4 || maxOrder > 1) {
      if (maxOrder === 3) return cls('notSp3', "This carbon is part of a triple bond, so it's linear (sp).");
      if (maxOrder === 2) return cls('notSp3', "This carbon is part of a double bond, so it's flat (sp2) and has only three groups.");
      return cls('notSp3', 'This carbon has only three groups.');
    }
    const nH = adj.filter(e => m.atoms[e.to].el === 'H').length;
    const hTot = Math.max(nH, a.hTotal || 0);
    if (hTot >= 2) {
      const what = hTot === 2 ? "two H's (it's a CH2)" : hTot === 3 ? "three H's (it's a CH3)" : "four H's (it's CH4)";
      return cls('multipleH', 'This carbon has ' + what + '. Two identical groups already rule it out.');
    }
    const ligs = adj.map(e => e.to);
    const groups = ligs.map(l => groupName(m, l, i));
    const r = rankSubstituents(m, i);
    if (r.capped) return cls('tooComplex', "This one is too tangled for me to rank. Try a smaller molecule.", { groups });
    if (!r.tie) {
      const list = joinList(r.order.map(l => groupName(m, l, i)));
      return cls('ok', 'Yes. {C} has four different groups: ' + list + '.', { groups });
    }
    // ties: pseudo / sameGroups / ringSame
    const busy = _busy || new Set();
    busy.add(i);
    const isRingPair = p => shareRing(m, i, p[0], p[1]) || (shareRing(m, i, p[0]) && shareRing(m, i, p[1]) &&
      subtree(m, i, p[0]).indexOf(p[1]) >= 0);
    for (const p of r.ties) {
      const pv = pseudoVerdict(m, i, p, coords, busy);
      if (pv === 'pseudo') {
        busy.delete(i);
        return cls('pseudo', "{C} sits between two groups that are built the same but differ in their own R/S. " +
          "That's a special case (often lowercase r/s) this tool doesn't cover.",
          { pair: p, comparison: compareLigands(m, i, p[0], p[1]), groups });
      }
      if (ezVerdict(m, i, p) === 'differ') {
        busy.delete(i);
        return cls('stereoDouble', "{C} has two groups that are built the same except for cis/trans (E/Z) at a double bond. " +
          "Z outranks E, so {C} is a real stereocenter, but that rule is a special case this tool doesn't cover.",
          { pair: p, comparison: compareLigands(m, i, p[0], p[1]), groups });
      }
    }
    busy.delete(i);
    const nonRing = r.ties.filter(p => !isRingPair(p));
    if (nonRing.length) {
      const p = nonRing[0];
      const A = groupName(m, p[0], i);
      return cls('sameGroups', 'This carbon has two identical groups (' + A + ' and ' + A + '). A stereocenter needs four different groups.',
        { pair: p, comparison: compareLigands(m, i, p[0], p[1]), groups });
    }
    const p = r.ties[0];
    let msg = shareRing(m, i, p[0], p[1])
      ? 'Going around the ring from this carbon both ways gives the same path, so two of its groups are identical.'
      : 'Two of the ring paths leaving this carbon are built the same, so two of its groups are identical.';
    if (ringHasOtherSubstituent(m, i, p)) msg += ' Cis/trans still matters on this ring, but it is named cis or trans, not R/S.';
    return cls('ringSame', msg, { pair: p, comparison: compareLigands(m, i, p[0], p[1]), groups });
  }

  function ringHasOtherSubstituent(m, i, pair) {
    const rings = (m.rings || []).filter(r => r.indexOf(i) >= 0 && r.indexOf(pair[0]) >= 0 && r.indexOf(pair[1]) >= 0);
    const ring = rings[0] || (m.rings || []).find(r => r.indexOf(i) >= 0);
    if (!ring) return false;
    return ring.some(j => j !== i && m.atoms[j].el === 'C' &&
      nbrs(m, j).some(t => m.atoms[t].el !== 'H' && ring.indexOf(t) < 0));
  }

  // For a tied pair: 'pseudo' if both branches carry stereocenters that are not identical, else null.
  function pseudoVerdict(m, i, pair, coords, busy) {
    const heavy = nHeavyOf(m);
    const centersIn = lig => {
      const sub0 = subtree(m, i, lig);
      const dist = bfsDist(m, i, lig);
      const out = [];
      for (const j of sub0) {
        if (j >= heavy || busy.has(j) || m.atoms[j].el !== 'C') continue;
        const c = classifyAtom(m, j, coords, busy);
        if (c.code === 'ok') out.push({ atom: j, d: dist[j] });
      }
      return out;
    };
    const ca = centersIn(pair[0]), cb = centersIn(pair[1]);
    if (!ca.length || !cb.length) return null;
    const lab = x => {
      const l = coords ? chirality(m, x.atom, coords) : chiralityFromTags(m, x.atom);
      return l;
    };
    const la = ca.map(x => ({ d: x.d, l: lab(x) })), lb = cb.map(x => ({ d: x.d, l: lab(x) }));
    if (la.some(x => !x.l) || lb.some(x => !x.l)) return 'pseudo';       // conservative without stereo info
    const sig = (arr, inv) => arr.map(x => x.d + ':' + (inv ? (x.l === 'R' ? 'S' : 'R') : x.l)).sort().join('|');
    if (sig(la, false) === sig(lb, false)) return null;                    // identical branches
    return 'pseudo';                                                       // mirror-image (or diastereomorphic) branches
  }
  // For a tied pair: 'differ' if the two branches differ only in the E/Z geometry of their double bonds
  // (CIP Rule 3, out of scope), else null. Uses bond.stereo = { x, y, cis } from NNChem's SMILES parser.
  function ezVerdict(m, i, pair) {
    const sig = lig => {
      const inSub = new Set(subtree(m, i, lig)), dist = bfsDist(m, i, lig), out = [];
      m.bonds.forEach(b => {
        if (b.order !== 2 || !b.stereo || !inSub.has(b.a) || !inSub.has(b.b)) return;
        const top = (end, other) => {
          const ns = nbrs(m, end).filter(t => t !== other);
          if (!ns.length) return null;
          if (ns.length === 1) return ns[0];
          const c = compareLigands(m, end, ns[0], ns[1]);
          return c.sign > 0 ? ns[0] : c.sign < 0 ? ns[1] : null;
        };
        const ta = top(b.a, b.b), tb = top(b.b, b.a);
        if (ta == null || tb == null) return;
        const xa = nbrs(m, b.a).indexOf(b.stereo.x) >= 0 ? b.stereo.x : b.stereo.y;
        const yb = xa === b.stereo.x ? b.stereo.y : b.stereo.x;
        let z = !!b.stereo.cis;
        if ((ta === xa) !== (tb === yb)) z = !z;
        out.push(Math.min(dist[b.a], dist[b.b]) + (z ? 'Z' : 'E'));
      });
      return out.sort().join('|');
    };
    const a = sig(pair[0]), b = sig(pair[1]);
    return (a || b) && a !== b ? 'differ' : null;
  }
  function bfsDist(m, block, start) {
    const adj = adjOf(m), d = {}; d[start] = 0;
    const q = [start], seen = new Set([block, start]);
    while (q.length) {
      const x = q.shift();
      for (const e of adj[x]) if (!seen.has(e.to)) { seen.add(e.to); d[e.to] = d[x] + 1; q.push(e.to); }
    }
    return d;
  }

  function findStereocenters(mol, coords) {
    const m = withH(mol);
    const heavy = nHeavyOf(m);
    const out = [];
    for (let i = 0; i < heavy; i++) {
      if (m.atoms[i].el !== 'C') continue;
      if (classifyAtom(m, i, coords).code !== 'ok') continue;
      out.push({ atom: i, label: chirality(m, i, coords), given: stereoGiven(m, i), ranking: rankSubstituents(m, i) });
    }
    return out;
  }

  /* ------------------------------------------------------------------ */
  /* Trace geometry (3.4)                                                 */
  /* ------------------------------------------------------------------ */
  function orientation(coords, center, fourAtom, toViewer) {
    const b4 = unit(sub(coords[fourAtom], coords[center]));
    const cos = dot(b4, unit(toViewer));
    return { four: cos <= -0.6 ? 'back' : cos >= 0.6 ? 'front' : 'side', cos };
  }
  function turn2D(c, p1, p2, p3) {
    const s = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
    if (Math.abs(s) < 1e-9) return null;
    return s > 0 ? 'cw' : 'ccw';
  }
  function labelFromTrace(turn, four, swaps) {
    if (turn !== 'cw' && turn !== 'ccw') return null;
    let r = turn === 'cw';
    if (four === 'front') r = !r;
    if ((swaps | 0) % 2) r = !r;
    return r ? 'R' : 'S';
  }
  function mirror(coords) { return coords.map(p => [-p[0], p[1], p[2]]); }

  /* ------------------------------------------------------------------ */
  /* Diagnosis (3.5)                                                      */
  /* ------------------------------------------------------------------ */
  function diagnose(mol, center, picked, expected, slot) {
    const m = withH(mol);
    const X = picked, Y = expected;
    const c = compareLigands(m, center, Y, X);
    const out = (code, msg) => ({ code, msg, comparison: c });
    if (m.atoms[X].el === 'H') return out('hFirst', 'H has the lowest atomic number (1), so it always comes last here.');
    const d = c.deciding;
    if (!d) return out('generic', '{Y} outranks {X} here, though the atomic-number rules alone do not separate them.');
    const a = d.a.el, b = d.b.el, Za = zOf(a), Zb = zOf(b);
    const setY = setText(c.boxes.a), setX = setText(c.boxes.b);
    if (b !== '0' && Za > Zb && EN[b] != null && EN[a] != null && EN[b] > EN[a]) {
      let msg = 'Rank by atomic number, not electronegativity: ' + a + ' (atomic number ' + Za + ') beats ' + b + ' (' + Zb + ').';
      if (c.sphere > 1) msg += ' Compare the highest atoms first, so one ' + a + ' beats any number of ' + b + "'s.";
      return out('electronegativity', msg);
    }
    if (c.sphere === 1) {
      return out('firstAtom', 'Look only at the atom attached to {C} first: ' + a + ' (' + Za + ') beats ' + b + ' (' + Zb +
        "). The rest of the group doesn't matter yet.");
    }
    const nd = compareLigands(m, center, Y, X, { duplicates: false });
    if (nd.sign <= 0) {
      return out('noDuplicate', 'Count a double bond as two bonds to the same atom (a triple bond as three). Then {Y} has ' +
        setY + ', which beats ' + setX + '.');
    }
    const heavyCount = lig => subtree(m, center, lig).filter(j => m.atoms[j].el !== 'H').length;
    if (heavyCount(X) > heavyCount(Y)) {
      return out('size', "Size doesn't decide it. Go out one atom at a time and stop at the first difference: {Y} has " +
        setY + ', {X} has ' + setX + ', and ' + a + ' beats ' + (b === '0' ? 'nothing' : b) + '.');
    }
    const sum = box => box.atoms.reduce((s, e) => s + zOf(e.el), 0);
    if (sum(c.boxes.b) > sum(c.boxes.a)) {
      return out('sum', "Don't add the atoms up. Compare highest to highest: " + a + ' beats ' + (b === '0' ? 'nothing' : b) +
        ', so ' + setY + ' wins over ' + setX + '.');
    }
    if (c.sphere >= 3) {
      return out('deeper', 'These two tie close to {C}. Move one atom further out, following the higher branch first: ' +
        setY + ' beats ' + setX + '.');
    }
    return out('generic', '{Y} outranks {X}: ' + setY + ' vs ' + setX + '. The first difference is ' + a + ' vs ' +
      (b === '0' ? 'nothing (0)' : b) + '.');
  }

  /* ------------------------------------------------------------------ */
  /* Meso and the summary (3.6)                                           */
  /* ------------------------------------------------------------------ */
  function normLabels(labels) {
    const out = {};
    if (Array.isArray(labels)) labels.forEach(x => { if (x && (x.label === 'R' || x.label === 'S')) out[x.atom] = x.label; });
    else Object.keys(labels || {}).forEach(k => { const v = labels[k]; if (v === 'R' || v === 'S') out[+k] = v; });
    return out;
  }

  function isMeso(mol, labels) {
    const L = normLabels(labels);
    const centers = Object.keys(L).map(Number);
    if (centers.length < 2) return false;
    const m = withH(mol);
    const n = nHeavyOf(m);
    const adj = adjOf(m);
    const hadj = [];
    for (let i = 0; i < n; i++) hadj.push(adj[i].filter(e => e.to < n));
    const hcount = i => adj[i].filter(e => m.atoms[e.to].el === 'H').length + (m.explicitH ? 0 : (m.atoms[i].hCount | 0));
    // invariant refinement
    let cl = [];
    for (let i = 0; i < n; i++) {
      const a = m.atoms[i];
      cl.push([a.el, a.charge | 0, hcount(i), hadj[i].length, hadj[i].map(e => e.order).sort().join('')].join('/'));
    }
    const compress = arr => { const u = Array.from(new Set(arr)).sort(); return arr.map(x => u.indexOf(x)); };
    let cls0 = compress(cl);
    for (let it = 0; it < n; it++) {
      const next = compress(cls0.map((c, i) => c + ':' + hadj[i].map(e => e.order + '-' + cls0[e.to]).sort().join(',')));
      const k1 = new Set(cls0).size, k2 = new Set(next).size;
      cls0 = next;
      if (k2 === k1) break;
    }
    const inv = l => (l === 'R' ? 'S' : 'R');
    const bond = [];
    for (let i = 0; i < n; i++) { bond.push({}); hadj[i].forEach(e => { bond[i][e.to] = e.order; }); }
    const sigma = new Array(n).fill(-1), used = new Array(n).fill(false);
    // order atoms: BFS so neighbors constrain early
    const orderA = [];
    const seen = new Array(n).fill(false);
    for (let s = 0; s < n; s++) {
      if (seen[s]) continue;
      seen[s] = true; const q = [s];
      while (q.length) { const x = q.shift(); orderA.push(x); hadj[x].forEach(e => { if (!seen[e.to]) { seen[e.to] = true; q.push(e.to); } }); }
    }
    let steps = 0;
    function ok(i, j) {
      if (cls0[i] !== cls0[j]) return false;
      if (L[i] != null) { if (L[j] == null || L[j] !== inv(L[i])) return false; }
      else if (L[j] != null) return false;
      for (const e of hadj[i]) {
        const sj = sigma[e.to];
        if (sj >= 0 && bond[j][sj] !== e.order) return false;
      }
      return true;
    }
    function bt(k) {
      if (++steps > 200000) return false;
      if (k === orderA.length) return true;
      const i = orderA[k];
      for (let j = 0; j < n; j++) {
        if (used[j] || !ok(i, j)) continue;
        sigma[i] = j; used[j] = true;
        if (bt(k + 1)) return true;
        sigma[i] = -1; used[j] = false;
      }
      return false;
    }
    return bt(0);
  }

  function describe(mol, labels, opts) {
    opts = opts || {};
    const L = normLabels(labels);
    const atoms = Object.keys(L).map(Number);
    const loc = opts.locants || null;
    const haveLoc = !!loc && atoms.every(a => loc[a] != null);
    atoms.sort((x, y) => (haveLoc ? (loc[x] - loc[y]) : 0) || (x - y));
    const list = atoms.map(a => ({ atom: a, label: L[a] }));
    const meso = atoms.length >= 2 && isMeso(mol, L);
    let descriptor = '', text = opts.base || '';
    if (atoms.length === 1) {
      descriptor = '(' + L[atoms[0]] + ')';
    } else if (atoms.length > 1 && haveLoc) {
      descriptor = '(' + atoms.map(a => loc[a] + L[a]).join(',') + ')';
    } else if (atoms.length > 1) {
      let lab;
      descriptor = atoms.map(a => {
        try { lab = C().atomLabel(mol, a); } catch (e) { lab = 'C' + (a + 1); }
        return lab + ' ' + L[a];
      }).join(', ');
    }
    if (opts.base && descriptor) {
      text = (atoms.length > 1 && !haveLoc) ? opts.base + ' (' + descriptor + ')' : descriptor + '-' + opts.base;
    } else if (!opts.base) text = descriptor;
    return { labels: list, meso, chiral: atoms.length > 0 && !meso, descriptor, text };
  }

  /* ------------------------------------------------------------------ */
  /* Plain-language walk through a comparison (extra, for the page)       */
  /* ------------------------------------------------------------------ */
  function explain(cmp, names) {
    names = names || {};
    const A = names.a || 'the first group', B = names.b || 'the second group';
    const out = [];
    if (!cmp || !cmp.rows) return out;
    cmp.rows.forEach(row => {
      const showA = row.a.filter((bx, k) => !bx.terminal || (row.diff && row.diff.box === k));
      const showB = row.b.filter((bx, k) => !bx.terminal || (row.diff && row.diff.box === k));
      if (row.sphere === 1) {
        const ea = row.a[0].atoms[0].el, eb = row.b[0].atoms[0].el;
        out.push(row.diff ? 'First atoms: ' + ea + ' vs ' + eb + '. ' + ea + ' has the higher atomic number, so ' + A + ' wins.'
          : 'First atoms: ' + ea + ' vs ' + eb + '. A tie, so look at what is attached.');
        return;
      }
      const ta = showA.map(setText).join(' '), tb = showB.map(setText).join(' ');
      if (!row.diff) {
        out.push('Next atoms out: ' + ta + ' vs ' + tb + '. Still a tie, so go one atom further.');
      } else {
        const ba = row.a[row.diff.box], bb = row.b[row.diff.box];
        const ea = ba.atoms[row.diff.pos].el, eb = bb.atoms[row.diff.pos].el;
        out.push('Next atoms out: ' + setText(ba) + ' vs ' + setText(bb) + '. First difference: ' + ea + ' beats ' +
          (eb === '0' ? 'nothing' : eb) + ', so ' + (cmp.sign > 0 ? A : B) + ' wins.');
      }
    });
    if (cmp.sign === 0 && cmp.capped) out.push('I stopped here; this comparison is too long to finish.');
    else if (cmp.sign === 0) out.push('No difference anywhere: the two groups are identical.');
    return out;
  }

  const api = {
    Z, EN,
    findStereocenters, classifyAtom, rankSubstituents, compareLigands,
    chirality, chiralityFromTags, stereoGiven,
    orientation, turn2D, labelFromTrace,
    diagnose, isMeso, mirror, describe, setText,
    explain
  };
  root.RSCip = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

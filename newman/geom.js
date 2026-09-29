/*
 * newman/geom.js  (window.NNGeom)
 * Newman Navigator: 3D embedding, force-field relaxation, bond rotation,
 * dihedral helpers, conformation naming and the torsional-energy teaching model.
 *
 * Contract: newman/SPEC.md section 4. Coordinates are in Angstrom, angles in degrees.
 * No dependencies besides NNChem (reached lazily inside functions). No Math.random:
 * the only randomness is a seeded mulberry32 (seed 20260929), so the same input gives
 * bit-identical coordinates on every run.
 *
 * Deviations / resolved ambiguities (documented here as the spec asks):
 *  1. Angle term: sp3 centres use theta0 = 112.0 deg when BOTH ends are heavy atoms
 *     (C-C-C, C-C-O ...) and 109.47 deg otherwise. With a pure 109.47 target the
 *     cyclohexane chair relaxes to ring torsions of ~60 deg and C-C-C angles of ~109.5.
 *     With 112 the chair comes out at ~57 deg (textbook ~55) and butane C-C-C at ~111.6.
 *  2. Linear (sp, 180 deg) angles use E = 2k(1 + cos theta), which matches k(theta - pi)^2
 *     near 180 deg but has no 1/sin(theta) singularity in its gradient.
 *  3. The sp2 out-of-plane term uses the signed volume V = (p1-c).((p2-c)x(p3-c)) instead of
 *     the exact height h. V is proportional to h near planarity; the force constant is
 *     scaled at setup (from the ideal bond lengths) so that E ~= 400 h^2 kJ/mol.
 *  4. Minimizer: L-BFGS (memory 8) with Armijo backtracking, a variant of the "gradient
 *     descent with backtracking" allowed by the spec, with a 0.3 A per-step cap. Same stop
 *     rules (max |grad| < 0.02, dE < 1e-7, maxIter 800).
 *  5. Naming: the spec's "ethane-like" case (staggered / eclipsed instead of anti /
 *     gauche) is applied whenever EITHER end is a symmetric rotor (all of its
 *     substituents, at least two of them, share one energy class, e.g. CH3, CCl3,
 *     C(CH3)3). Every staggered position is then equivalent, so "anti" or
 *     "gauche" would be misleading (propane C1-C2, ethanol C1-C2, butane C1-C2).
 *  6. Profile.terms entries carry {front, back, ecl, gauche} with front/back = energy
 *     class ('H','Me','iPr','tBu','F','Cl','Br','I','O','N','S'), plus a display
 *     `label` like 'CH3/CH3' and `oh` (true when an O-H hydrogen-bond variant was used).
 *  7. The profile zero is the minimum of the 2-degree grid (so min(points.kJ) is exactly
 *     0); energyAt uses the same shift and is clamped at 0 for angles between samples.
 *  8. ConfName for ring bonds whose ring is neither chair, boat nor twist-boat: key
 *     'ring', text 'Ring bond (half-chair)' etc. (short form before '(' = 'Ring bond').
 *     Triple bonds and sp ends return key 'fixed' with their own text.
 *  9. If NNChem lacks a helper (hybrid, groupInfo, EL, rings), a small internal fallback
 *     is used so this module never throws on partial input.
 */
(function (root) {
  'use strict';

  // ---------------------------------------------------------------- deps (lazy)
  function chem() {
    if (root.NNChem) return root.NNChem;
    if (typeof require === 'function') {
      try { return require('./chem.js'); } catch (e) { /* not available */ }
    }
    return null;
  }

  // ---------------------------------------------------------------- constants
  const DEG = Math.PI / 180;
  const SEED = 20260929;
  const TETRA = 109.4712206;       // ideal sp3 angle
  const TETRA_HEAVY = 112.0;       // sp3 angle between two heavy substituents (deviation 1)
  const KB = 1400, KTH = 250, KIMP = 400, KREP = 25, SIGMA_F = 0.82;
  const V3_HH = 1.2, V3_HEAVY = 1.8, K_PLANAR = 20;

  const VDW_FALLBACK = { H: 1.20, C: 1.70, N: 1.55, O: 1.52, F: 1.47, Cl: 1.75, Br: 1.85, I: 1.98, S: 1.80, B: 1.92, P: 1.80 };
  const COV = { H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57, Cl: 1.02, Br: 1.20, I: 1.39, S: 1.05, B: 0.84, P: 1.07 };
  // Bond lengths (A), keyed by sorted element pair + order.
  const R0 = {
    'C-C-1': 1.53, 'C-C-2': 1.34, 'C-C-3': 1.20, 'C-H-1': 1.09,
    'C-N-1': 1.47, 'C-N-2': 1.28, 'C-N-3': 1.16,
    'C-O-1': 1.43, 'C-O-2': 1.21, 'H-O-1': 0.96, 'H-N-1': 1.01, 'H-S-1': 1.34,
    'C-F-1': 1.35, 'C-Cl-1': 1.77, 'Br-C-1': 1.94, 'C-I-1': 2.14, 'C-S-1': 1.82,
    'N-N-1': 1.45, 'N-O-1': 1.40, 'O-O-1': 1.48, 'B-C-1': 1.56, 'C-P-1': 1.84
  };

  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // ---------------------------------------------------------------- vector helpers
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const scl = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = (a) => Math.sqrt(dot(a, a));
  function unit(a) { const l = len(a); return l > 1e-12 ? scl(a, 1 / l) : [0, 0, 0]; }
  function perp(v, axis) { return sub(v, scl(axis, dot(v, axis))); }       // axis must be unit
  function anyPerp(a) {                                                     // some unit vector perpendicular to unit a
    const t = Math.abs(a[2]) < 0.9 ? [0, 0, 1] : [0, 1, 0];
    return unit(perp(t, a));
  }
  // Rodrigues rotation of v about unit axis k by angle (rad)
  function rotVec(v, k, ang) {
    const c = Math.cos(ang), s = Math.sin(ang);
    const kxv = cross(k, v), kdv = dot(k, v) * (1 - c);
    return [v[0] * c + kxv[0] * s + k[0] * kdv, v[1] * c + kxv[1] * s + k[1] * kdv, v[2] * c + kxv[2] * s + k[2] * kdv];
  }
  // 3x3 rotation matrix (row-major) about unit axis k by ang
  function rotMat(k, ang) {
    const c = Math.cos(ang), s = Math.sin(ang), t = 1 - c, [x, y, z] = k;
    return [
      [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
      [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
      [t * x * z - s * y, t * y * z + s * x, t * z * z + c]
    ];
  }
  const mulMV = (m, v) => [dot(m[0], v), dot(m[1], v), dot(m[2], v)];
  function matMul(a, b) {
    const r = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    return r;
  }
  // Rotation taking unit vector a onto unit vector b
  function alignMat(a, b) {
    const c = dot(a, b);
    if (c > 1 - 1e-12) return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    if (c < -1 + 1e-12) return rotMat(anyPerp(a), Math.PI);
    return rotMat(unit(cross(a, b)), Math.acos(c));
  }
  function wrap180(d) { d = ((d % 360) + 360) % 360; return d > 180 ? d - 360 : d; }      // (-180, 180]
  function wrap360(d) { d = ((d % 360) + 360) % 360; return d >= 360 ? 0 : d; }            // [0, 360)
  function det3(a, b, c) { return dot(a, cross(b, c)); }

  // ---------------------------------------------------------------- graph helpers
  function adjacency(mol) {
    const n = mol.atoms.length, adj = [];
    for (let i = 0; i < n; i++) adj.push([]);
    mol.bonds.forEach((b, bi) => {
      adj[b.a].push({ n: b.b, order: b.order || 1, bond: bi });
      adj[b.b].push({ n: b.a, order: b.order || 1, bond: bi });
    });
    adj.forEach(l => l.sort((p, q) => p.n - q.n));
    return adj;
  }
  function nbrs(adj, i) { return adj[i].map(e => e.n); }
  function isH(mol, i) { return mol.atoms[i].el === 'H'; }

  function hybridOf(mol, adj, i) {
    const C = chem();
    if (C && typeof C.hybrid === 'function') {
      try { const h = C.hybrid(mol, i); if (h === 'sp' || h === 'sp2' || h === 'sp3') return h; } catch (e) { /* fall through */ }
    }
    let dbl = 0, tpl = 0;
    adj[i].forEach(e => { if (e.order === 2) dbl++; else if (e.order === 3) tpl++; });
    if (tpl || dbl >= 2) return 'sp';
    if (dbl === 1) return 'sp2';
    return 'sp3';
  }

  function vdwOf(el) {
    const C = chem();
    const e = C && C.EL && C.EL[el];
    return (e && typeof e.vdw === 'number') ? e.vdw : (VDW_FALLBACK[el] || 1.7);
  }

  function bondR0(e1, e2, order) {
    const k = [e1, e2].sort().join('-') + '-' + order;
    if (R0[k] !== undefined) return R0[k];
    return (COV[e1] || 0.76) + (COV[e2] || 0.76) - 0.10 * (order - 1);
  }

  // Ring list (heavy-atom cycles). Prefer NNChem's perception; otherwise a small
  // BFS smallest-cycle-per-ring-bond fallback.
  function ringsOf(mol, adj) {
    if (Array.isArray(mol.rings)) return mol.rings;
    const rings = [], seen = new Set();
    mol.bonds.forEach(b => {
      if (isH(mol, b.a) || isH(mol, b.b)) return;
      // shortest path a -> b avoiding the direct bond
      const prev = new Map([[b.a, -1]]), q = [b.a];
      while (q.length) {
        const x = q.shift();
        if (x === b.b) break;
        for (const e of adj[x]) {
          if (isH(mol, e.n) || prev.has(e.n)) continue;
          if ((x === b.a && e.n === b.b)) continue;
          prev.set(e.n, x); q.push(e.n);
        }
      }
      if (!prev.has(b.b)) return;
      const cyc = [];
      for (let x = b.b; x !== -1; x = prev.get(x)) cyc.push(x);
      const key = cyc.slice().sort((p, q2) => p - q2).join(',');
      if (!seen.has(key)) { seen.add(key); rings.push(cyc.reverse()); }
    });
    return rings;
  }

  // Is bond a-b on a cycle? (b reachable from a without using the bond)
  function bondInCycle(adj, a, b) {
    const seen = new Set([a]), q = [a];
    while (q.length) {
      const x = q.shift();
      for (const e of adj[x]) {
        if (x === a && e.n === b) continue;
        if (e.n === b) return true;
        if (!seen.has(e.n)) { seen.add(e.n); q.push(e.n); }
      }
    }
    return false;
  }

  // Group info with a tiny fallback (cls/score only) when NNChem.groupInfo is missing.
  function groupOf(mol, atomIdx, fromIdx) {
    const C = chem();
    if (C && typeof C.groupInfo === 'function') {
      try { const g = C.groupInfo(mol, atomIdx, fromIdx); if (g) return g; } catch (e) { /* fall through */ }
    }
    const a = mol.atoms[atomIdx], el = a.el;
    if (el === 'H') return { label: 'H', cls: 'H', score: 0, el: 'H', ring: false, hasH: false };
    const adj = adjacency(mol);
    const heavyN = adj[atomIdx].filter(e => e.n !== fromIdx && !isH(mol, e.n)).length;
    const hasH = adj[atomIdx].some(e => isH(mol, e.n)) || (a.hCount || 0) > 0;
    const hetero = { F: 2.0, Cl: 2.6, Br: 2.8, I: 2.9, O: 2.2, N: 2.3, S: 2.7 };
    if (hetero[el] !== undefined) return { label: el, cls: el, score: hetero[el], el, ring: false, hasH };
    const cls = heavyN >= 3 ? 'tBu' : heavyN === 2 ? 'iPr' : 'Me';
    const score = cls === 'tBu' ? 5 : cls === 'iPr' ? 4 : 3 + 0.01 * (1 + heavyN);
    return { label: el, cls, score, el, ring: false, hasH };
  }

  // ================================================================ build3D
  function build3D(mol, opts) {
    opts = opts || {};
    const C = chem();
    let m;
    if (mol.explicitH) m = JSON.parse(JSON.stringify(mol));
    else if (C && typeof C.addHydrogens === 'function') m = C.addHydrogens(mol);
    else throw new Error('NNGeom.build3D needs NNChem.addHydrogens');
    const rand = mulberry32(opts.seed === undefined ? SEED : opts.seed);
    const coords = embed(m, opts);
    for (let i = 0; i < coords.length; i++) {
      for (let d = 0; d < 3; d++) coords[i][d] += (rand() - 0.5) * 0.02;       // +-0.01 A jitter
    }
    return { mol: m, coords: relax(m, coords, { maxIter: opts.maxIter }) };
  }

  // Initial geometry: ring templates + BFS placement with ideal directions (SPEC 4.3)
  function embed(m, eopts) {
    eopts = eopts || {};
    const n = m.atoms.length, adj = adjacency(m), rings = ringsOf(m, adj);
    const hyb = []; for (let i = 0; i < n; i++) hyb.push(hybridOf(m, adj, i));
    const X = new Array(n).fill(null);
    const queue = [];
    const ringPlaced = rings.map(() => false);
    const ringsOfAtom = []; for (let i = 0; i < n; i++) ringsOfAtom.push([]);
    rings.forEach((r, ri) => r.forEach(a => ringsOfAtom[a].push(ri)));
    const r0 = (i, j) => { const e = adj[i].find(q => q.n === j); return bondR0(m.atoms[i].el, m.atoms[j].el, e ? e.order : 1); };

    // Score used to choose the reference ("biggest") substituent for all-anti placement
    const scoreCache = new Map();
    function score(x, from) {
      const k = x * 100003 + from;
      if (!scoreCache.has(k)) scoreCache.set(k, groupOf(m, x, from).score || 0);
      return scoreCache.get(k);
    }

    // Local template coordinates for a ring (centred at origin, mean plane z = 0)
    function ringTemplate(ri) {
      const ring = rings[ri], k = ring.length;
      const allSp3 = ring.every(a => hyb[a] === 'sp3');
      let L = 0;
      for (let j = 0; j < k; j++) L += r0(ring[j], ring[(j + 1) % k]);
      L /= k;
      if (k === 6 && allSp3) {
        // Chair: adjacent atoms 60 deg apart at radius rho, heights +-h, with C-C-C 111 deg.
        const th = 111 * DEG;
        const rho = Math.sqrt((2 * L * L * (1 - Math.cos(th))) / 3);
        const h = Math.sqrt(Math.max(0, L * L - rho * rho)) / 2;
        return ring.map((_, j) => [rho * Math.cos(j * Math.PI / 3), rho * Math.sin(j * Math.PI / 3), j % 2 ? -h : h]);
      }
      const rho = L / (2 * Math.sin(Math.PI / k));
      const allSp2 = ring.every(a => hyb[a] !== 'sp3');
      const puck = (k <= 3 || allSp2) ? 0 : k === 4 ? 0.15 : 0.25;
      return ring.map((_, j) => [rho * Math.cos(2 * Math.PI * j / k), rho * Math.sin(2 * Math.PI * j / k), j % 2 ? -puck : puck]);
    }

    // Exocyclic direction for ring atom j of a template (equatorial for sp3)
    function templateExo(T, j, sp3) {
      const k = T.length, p = T[j];
      const b1 = unit(sub(T[(j + k - 1) % k], p)), b2 = unit(sub(T[(j + 1) % k], p));
      const d = unit(scl(add(b1, b2), -1));
      if (!sp3) return d;
      const nn = unit(cross(b1, b2)), c = Math.cos(54.7356 * DEG), s = Math.sin(54.7356 * DEG);
      const e1 = add(scl(d, c), scl(nn, s)), e2 = sub(scl(d, c), scl(nn, s));
      return Math.abs(e1[2]) < Math.abs(e2[2]) ? e1 : e2;   // equatorial = closest to mean plane
    }

    // Place ring ri as a template: either free-standing (origin) or hanging off atom c,
    // already placed and bonded to placed atom P (template exo bond of c points at P).
    function placeRing(ri, c, P) {
      const ring = rings[ri], T = ringTemplate(ri);
      ringPlaced[ri] = true;
      if (c === undefined) {
        ring.forEach((a, j) => { X[a] = T[j].slice(); });
        ring.forEach(a => queue.push(a));
        return;
      }
      const j = ring.indexOf(c);
      const exo = templateExo(T, j, hyb[c] === 'sp3');
      const target = unit(sub(X[P], X[c]));
      let R = alignMat(exo, target);
      // Free spin about the c->P axis: put the ring centroid anti to P's biggest other substituent.
      const cen = mulMV(R, scl(T[j], -1));   // centroid (origin) relative to c
      const others = adj[P].map(e => e.n).filter(x => x !== c && X[x]);
      if (others.length) {
        let best = others[0];
        others.forEach(x => { if (score(x, P) > score(best, P) + 1e-9) best = x; });
        const want = unit(scl(perp(sub(X[best], X[P]), target), -1));
        const have = unit(perp(cen, target));
        if (len(want) > 0 && len(have) > 0) {
          const ang = Math.atan2(dot(cross(have, want), target), dot(have, want));
          R = matMul(rotMat(target, ang), R);
        }
      }
      ring.forEach((a, jj) => {
        if (a === c) return;
        X[a] = add(X[c], mulMV(R, sub(T[jj], T[j])));
      });
      ring.forEach(a => { if (a !== c) queue.push(a); });
    }

    function ringNormal(ri) {
      const ring = rings[ri]; let nn = [0, 0, 0];
      const cen = ring.reduce((s, a) => add(s, X[a] || [0, 0, 0]), [0, 0, 0]).map(v => v / ring.length);
      for (let j = 0; j < ring.length; j++) {
        const p = X[ring[j]], q = X[ring[(j + 1) % ring.length]];
        if (p && q) nn = add(nn, cross(sub(p, cen), sub(q, cen)));
      }
      return unit(nn);
    }

    // Place all unplaced neighbours of P.
    function placeChildren(P) {
      const all = adj[P].map(e => e.n);
      let un = all.filter(x => !X[x]);
      if (!un.length) return;
      // heavy children first, biggest group first (it takes the anti slot, psi = 0), then H
      un.sort((a, b) => (isH(m, a) - isH(m, b)) || (score(b, P) - score(a, P)) || (a - b));
      let placed = all.filter(x => X[x]);
      if (!placed.length) {
        // root atom: first child along +x, then treat as the one-placed case
        const c0 = un[0];
        X[c0] = add(X[P], [r0(P, c0), 0, 0]);
        afterPlace(c0, P);
        un = un.slice(1);
        placed = [c0];
        if (!un.length) return;
      }
      const bvec = placed.map(x => unit(sub(X[x], X[P])));
      let dirs = [];
      const h = hyb[P];
      if (placed.length === 1) {
        const G = placed[0], a = bvec[0];
        // reference: biggest other substituent on G (for the all-anti rule)
        const gOthers = adj[G].map(e => e.n).filter(x => x !== P && X[x]);
        // cis/trans double bond G=P with a stored E/Z: the reference is the marked
        // substituent on G, and P's marked substituent goes trans (psi 0) or cis (psi 180).
        const gpBond = m.bonds[bondIndex(m, G, P)];
        const ez = gpBond && gpBond.order === 2 && gpBond.stereo ? gpBond.stereo : null;
        let ezRef = -1, ezMine = -1;
        if (ez) {
          const onG = [ez.x, ez.y].find(q => q !== P && adj[G].some(e => e.n === q));
          const onP = [ez.x, ez.y].find(q => q !== G && adj[P].some(e => e.n === q));
          if (onG != null && X[onG] && onP != null && un.includes(onP)) {
            ezRef = onG; ezMine = onP;
            const other = un.filter(q => q !== onP);
            un = ez.cis ? [other[0], onP].concat(other.slice(1)).filter(q => q != null) : [onP].concat(other);
          }
        }
        let eX;
        if (gOthers.length) {
          let best = gOthers[0];
          gOthers.forEach(x => { if (score(x, G) > score(best, G) + 1e-9) best = x; });
          if (ezRef >= 0) best = ezRef;
          eX = unit(perp(sub(X[best], X[G]), a));
          if (len(eX) < 0.5) eX = anyPerp(a);
        } else eX = anyPerp(a);
        const e1 = scl(eX, -1), e2 = cross(a, e1);
        let ang, psis;
        if (h === 'sp') { ang = 180; psis = [0]; }
        else if (h === 'sp2') { ang = 120; psis = [0, 180]; }
        else { ang = TETRA; psis = [0, 120, 240]; }
        dirs = psis.map(ps => unit(add(scl(a, Math.cos(ang * DEG)), scl(add(scl(e1, Math.cos(ps * DEG)), scl(e2, Math.sin(ps * DEG))), Math.sin(ang * DEG)))));
        if (ezMine >= 0 && ez.cis && un.length === 1 && dirs.length > 1) dirs = [dirs[1], dirs[0]];
      } else {
        const sum = bvec.reduce(add, [0, 0, 0]);
        if (h === 'sp3' && placed.length === 2) {
          const d = unit(scl(sum, -1)), nn = unit(cross(bvec[0], bvec[1]));
          const c = Math.cos(54.7356 * DEG), s = Math.sin(54.7356 * DEG);
          dirs = [add(scl(d, c), scl(nn, s)), sub(scl(d, c), scl(nn, s))];
          // ring atom: heavy substituent goes equatorial (direction closest to the ring plane)
          const ri = ringsOfAtom[P].find(r => ringPlaced[r]);
          if (ri !== undefined) {
            const rn = ringNormal(ri);
            const eq1 = Math.abs(dot(dirs[1], rn)) < Math.abs(dot(dirs[0], rn));
            // eopts.axial (chair flip): the heavy substituent goes axial instead
            if (eq1 !== !!eopts.axial) dirs = [dirs[1], dirs[0]];
          }
        } else if (len(sum) > 1e-6) {
          dirs = [unit(scl(sum, -1))];
        }
      }
      // More children than ideal slots (odd valences): spread the extras deterministically.
      while (dirs.length < un.length) {
        const s = bvec.concat(dirs).reduce(add, [0, 0, 0]);
        let d = unit(scl(s, -1));
        if (len(d) < 0.5) d = anyPerp(bvec[0] || [1, 0, 0]);
        dirs.push(unit(add(d, scl(anyPerp(d), 0.3 * dirs.length))));
      }
      const assign = un.map((x, i) => dirs[i]);
      // Chirality: swap the last two children if the determinant sign is wrong (SPEC 4.3.5)
      const at = m.atoms[P];
      if (at.chiral && Array.isArray(at.chiralNbrs) && at.chiralNbrs.length === 4 && un.length >= 2 && at.chiralNbrs.every(q => q >= 0)) {
        const pos = (x) => { const k = un.indexOf(x); return k >= 0 ? add(X[P], assign[k]) : X[x]; };
        const p = at.chiralNbrs.map(pos);
        if (p.every(Boolean)) {
          const d = det3(sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
          const want = at.chiral === '@' ? -1 : 1;
          if (Math.sign(d) !== want) {
            const k = un.length;
            const t = assign[k - 1]; assign[k - 1] = assign[k - 2]; assign[k - 2] = t;
          }
        }
      }
      un.forEach((x, i) => {
        X[x] = add(X[P], scl(assign[i], r0(P, x)));
        afterPlace(x, P);
      });
    }

    function afterPlace(x, P) {
      // If x starts a ring nobody has placed yet, drop the whole ring template in now.
      for (const ri of ringsOfAtom[x]) {
        if (ringPlaced[ri]) continue;
        if (rings[ri].every(a => a === x || !X[a])) { placeRing(ri, x, P); return; }
      }
      queue.push(x);
    }

    // Components (normally one). Start from the first ring, else heavy atom 0.
    let offset = 0;
    for (;;) {
      let seed = -1, seedRing = -1;
      for (let ri = 0; ri < rings.length; ri++) if (!ringPlaced[ri] && rings[ri].every(a => !X[a])) { seedRing = ri; break; }
      if (seedRing < 0) for (let i = 0; i < n; i++) if (!X[i]) { seed = i; break; }
      if (seedRing < 0 && seed < 0) break;
      if (seedRing >= 0) {
        placeRing(seedRing);
        if (offset) rings[seedRing].forEach(a => { X[a][0] += offset; });
      } else {
        X[seed] = [offset, 0, 0];
        queue.push(seed);
      }
      while (queue.length) placeChildren(queue.shift());
      offset += 8;
      if (offset > 8 * n + 8) break;   // safety
    }
    return X.map(p => p || [0, 0, 0]);
  }

  // ================================================================ force field
  // Precomputed term lists, cached per mol object.
  const ffCache = typeof WeakMap !== 'undefined' ? new WeakMap() : null;

  function setupFF(m) {
    if (ffCache && ffCache.has(m)) return ffCache.get(m);
    const n = m.atoms.length, adj = adjacency(m), rings = ringsOf(m, adj);
    const hyb = []; for (let i = 0; i < n; i++) hyb.push(hybridOf(m, adj, i));
    const el = m.atoms.map(a => a.el);
    const ff = { n, bonds: [], angles: [], tors: [], imps: [], nb: [] };

    m.bonds.forEach(b => ff.bonds.push([b.a, b.b, bondR0(el[b.a], el[b.b], b.order || 1)]));

    const smallRing = (c, i, k) => {        // size of the smallest 3/4-ring containing all three
      let best = 0;
      rings.forEach(r => { if (r.length <= 4 && r.includes(c) && r.includes(i) && r.includes(k)) best = best ? Math.min(best, r.length) : r.length; });
      return best;
    };
    for (let c = 0; c < n; c++) {
      const nb = nbrs(adj, c);
      for (let p = 0; p < nb.length; p++) for (let q = p + 1; q < nb.length; q++) {
        const i = nb[p], k = nb[q];
        let th0;
        if (hyb[c] === 'sp') th0 = 180;
        else if (hyb[c] === 'sp2') th0 = 120;
        else th0 = (el[i] !== 'H' && el[k] !== 'H') ? TETRA_HEAVY : TETRA;
        const sr = smallRing(c, i, k);
        if (sr === 3) th0 = 60; else if (sr === 4) th0 = 90;
        ff.angles.push([i, c, k, th0 * DEG, th0 >= 179.9]);
      }
    }

    // Torsions (3-fold, sp3-sp3 single bonds) and planarity (sp2=sp2 double bonds)
    m.bonds.forEach(b => {
      const A = b.a, B = b.b, order = b.order || 1;
      const XA = nbrs(adj, A).filter(x => x !== B), YB = nbrs(adj, B).filter(y => y !== A);
      if (order === 1 && hyb[A] === 'sp3' && hyb[B] === 'sp3') {
        XA.forEach(x => YB.forEach(y => {
          const V3 = (el[x] === 'H' && el[y] === 'H') ? V3_HH : V3_HEAVY;
          ff.tors.push([x, A, B, y, V3 / 2, 1, 3]);       // K(1 + cos 3phi)
        }));
      } else if (order === 1 && ((hyb[A] === 'sp3' && hyb[B] === 'sp2') || (hyb[A] === 'sp2' && hyb[B] === 'sp3'))) {
        // allylic / carbonyl single bond: X-C(sp3)-C(sp2)=Y prefers X eclipsing the double bond,
        // E = K(1 - cos 3phi) with K ~ 1.4 kJ/mol for every X on the sp3 end
        const s3 = hyb[A] === 'sp3' ? A : B, s2 = s3 === A ? B : A;
        const Xs = nbrs(adj, s3).filter(x => x !== s2);
        adj[s2].filter(e => e.n !== s3 && e.order === 2).forEach(e => {
          Xs.forEach(x => ff.tors.push([x, s3, s2, e.n, 1.4, -1, 3]));
        });
      } else if (order === 2 && hyb[A] === 'sp2' && hyb[B] === 'sp2') {
        XA.forEach(x => YB.forEach(y => ff.tors.push([x, A, B, y, K_PLANAR, -1, 2])));   // K(1 - cos 2phi)
      }
    });

    // sp2 out-of-plane (signed volume, scaled so E ~= KIMP * h^2; deviation 3)
    for (let c = 0; c < n; c++) {
      if (hyb[c] !== 'sp2') continue;
      const nb = nbrs(adj, c);
      if (nb.length !== 3) continue;
      const L = nb.reduce((s, x) => s + bondR0(el[c], el[x], (adj[c].find(e => e.n === x) || {}).order || 1), 0) / 3;
      const twoArea = 1.5 * Math.sqrt(3) * L * L;
      ff.imps.push([c, nb[0], nb[1], nb[2], KIMP / (twoArea * twoArea)]);
    }

    // Non-bonded repulsion for pairs >= 3 bonds apart (1-4 scaled 0.5)
    const topo = [];
    for (let i = 0; i < n; i++) {
      const d = new Map([[i, 0]]); let front = [i];
      for (let depth = 1; depth <= 3; depth++) {
        const next = [];
        front.forEach(x => adj[x].forEach(e => { if (!d.has(e.n)) { d.set(e.n, depth); next.push(e.n); } }));
        front = next;
      }
      topo.push(d);
    }
    const vdw = el.map(vdwOf);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const d = topo[i].has(j) ? topo[i].get(j) : 99;
      if (d < 3) continue;
      ff.nb.push([i, j, SIGMA_F * (vdw[i] + vdw[j]), d === 3 ? KREP * 0.5 : KREP]);
    }
    if (ffCache) ffCache.set(m, ff);
    return ff;
  }

  // Energy and gradient (g is overwritten). x is a flat Float64Array of 3n.
  function energyGrad(ff, x, g) {
    g.fill(0);
    let E = 0;
    // bonds
    for (const [i, j, r0] of ff.bonds) {
      const dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-9, dr = r - r0;
      E += KB * dr * dr;
      const f = 2 * KB * dr / r;
      g[3 * i] += f * dx; g[3 * i + 1] += f * dy; g[3 * i + 2] += f * dz;
      g[3 * j] -= f * dx; g[3 * j + 1] -= f * dy; g[3 * j + 2] -= f * dz;
    }
    // angles
    for (const [i, c, k, th0, linear] of ff.angles) {
      const ux = x[3 * i] - x[3 * c], uy = x[3 * i + 1] - x[3 * c + 1], uz = x[3 * i + 2] - x[3 * c + 2];
      const vx = x[3 * k] - x[3 * c], vy = x[3 * k + 1] - x[3 * c + 1], vz = x[3 * k + 2] - x[3 * c + 2];
      const lu = Math.sqrt(ux * ux + uy * uy + uz * uz) || 1e-9, lv = Math.sqrt(vx * vx + vy * vy + vz * vz) || 1e-9;
      let cs = (ux * vx + uy * vy + uz * vz) / (lu * lv);
      cs = Math.max(-1, Math.min(1, cs));
      let dEdc;
      if (linear) {                         // E = 2k(1 + cos)  (deviation 2)
        E += 2 * KTH * (1 + cs);
        dEdc = 2 * KTH;
      } else {
        const th = Math.acos(cs), dth = th - th0;
        E += KTH * dth * dth;
        const sn = Math.max(Math.sqrt(1 - cs * cs), 1e-8);
        dEdc = 2 * KTH * dth * (-1 / sn);
      }
      // dc/du = v/(lu lv) - c u/lu^2 ; dc/dv = u/(lu lv) - c v/lv^2
      const a1 = 1 / (lu * lv), bu = cs / (lu * lu), bv = cs / (lv * lv);
      const gix = dEdc * (vx * a1 - ux * bu), giy = dEdc * (vy * a1 - uy * bu), giz = dEdc * (vz * a1 - uz * bu);
      const gkx = dEdc * (ux * a1 - vx * bv), gky = dEdc * (uy * a1 - vy * bv), gkz = dEdc * (uz * a1 - vz * bv);
      g[3 * i] += gix; g[3 * i + 1] += giy; g[3 * i + 2] += giz;
      g[3 * k] += gkx; g[3 * k + 1] += gky; g[3 * k + 2] += gkz;
      g[3 * c] -= gix + gkx; g[3 * c + 1] -= giy + gky; g[3 * c + 2] -= giz + gkz;
    }
    // torsions / planarity: E = K (1 + s cos(n phi))
    for (const [i, j, k, l, K, s, nn] of ff.tors) {
      const Fx = x[3 * i] - x[3 * j], Fy = x[3 * i + 1] - x[3 * j + 1], Fz = x[3 * i + 2] - x[3 * j + 2];
      const Gx = x[3 * j] - x[3 * k], Gy = x[3 * j + 1] - x[3 * k + 1], Gz = x[3 * j + 2] - x[3 * k + 2];
      const Hx = x[3 * l] - x[3 * k], Hy = x[3 * l + 1] - x[3 * k + 1], Hz = x[3 * l + 2] - x[3 * k + 2];
      const Ax = Fy * Gz - Fz * Gy, Ay = Fz * Gx - Fx * Gz, Az = Fx * Gy - Fy * Gx;
      const Bx = Hy * Gz - Hz * Gy, By = Hz * Gx - Hx * Gz, Bz = Hx * Gy - Hy * Gx;
      const A2 = Ax * Ax + Ay * Ay + Az * Az, B2 = Bx * Bx + By * By + Bz * Bz;
      const Gl = Math.sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
      if (A2 < 1e-10 || B2 < 1e-10 || Gl < 1e-8) continue;
      const AB = Math.sqrt(A2 * B2);
      const cphi = (Ax * Bx + Ay * By + Az * Bz) / AB;
      // (B x A) . G
      const BxAx = By * Az - Bz * Ay, BxAy = Bz * Ax - Bx * Az, BxAz = Bx * Ay - By * Ax;
      const sphi = (BxAx * Gx + BxAy * Gy + BxAz * Gz) / (AB * Gl);
      const phi = Math.atan2(sphi, cphi);
      E += K * (1 + s * Math.cos(nn * phi));
      const dEdphi = -K * s * nn * Math.sin(nn * phi);
      const FG = Fx * Gx + Fy * Gy + Fz * Gz, HG = Hx * Gx + Hy * Gy + Hz * Gz;
      const fa = -Gl / A2, fb = Gl / B2;           // dphi/dri = fa*A ; dphi/drl = fb*B
      const ta = FG / (A2 * Gl), tb = HG / (B2 * Gl);
      for (let d = 0; d < 3; d++) {
        const Ad = d === 0 ? Ax : d === 1 ? Ay : Az, Bd = d === 0 ? Bx : d === 1 ? By : Bz;
        const gi = fa * Ad, gl = fb * Bd;
        const gj = -gi + ta * Ad - tb * Bd;
        const gk = -gl - ta * Ad + tb * Bd;
        g[3 * i + d] += dEdphi * gi; g[3 * l + d] += dEdphi * gl;
        g[3 * j + d] += dEdphi * gj; g[3 * k + d] += dEdphi * gk;
      }
    }
    // impropers: E = kv * V^2, V = (p1-c).((p2-c) x (p3-c))
    for (const [c, p1, p2, p3, kv] of ff.imps) {
      const a = [x[3 * p1] - x[3 * c], x[3 * p1 + 1] - x[3 * c + 1], x[3 * p1 + 2] - x[3 * c + 2]];
      const b = [x[3 * p2] - x[3 * c], x[3 * p2 + 1] - x[3 * c + 1], x[3 * p2 + 2] - x[3 * c + 2]];
      const d = [x[3 * p3] - x[3 * c], x[3 * p3 + 1] - x[3 * c + 1], x[3 * p3 + 2] - x[3 * c + 2]];
      const bxd = cross(b, d), dxa = cross(d, a), axb = cross(a, b);
      const V = dot(a, bxd);
      E += kv * V * V;
      const f = 2 * kv * V;
      for (let t = 0; t < 3; t++) {
        g[3 * p1 + t] += f * bxd[t]; g[3 * p2 + t] += f * dxa[t]; g[3 * p3 + t] += f * axb[t];
        g[3 * c + t] -= f * (bxd[t] + dxa[t] + axb[t]);
      }
    }
    // repulsion
    for (const [i, j, sg, k] of ff.nb) {
      const dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1], dz = x[3 * i + 2] - x[3 * j + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= sg * sg) continue;
      const r = Math.sqrt(r2) || 1e-9, dr = sg - r;
      E += k * dr * dr;
      const f = -2 * k * dr / r;
      g[3 * i] += f * dx; g[3 * i + 1] += f * dy; g[3 * i + 2] += f * dz;
      g[3 * j] -= f * dx; g[3 * j + 1] -= f * dy; g[3 * j + 2] -= f * dz;
    }
    return E;
  }

  // L-BFGS with Armijo backtracking (deviation 4)
  function relax(m, coords, opts) {
    opts = opts || {};
    const maxIter = opts.maxIter || 800;
    const ff = setupFF(m), n = ff.n, N = 3 * n;
    const frozen = new Uint8Array(n);
    (opts.frozen || []).forEach(i => { if (i >= 0 && i < n) frozen[i] = 1; });
    const x = new Float64Array(N);
    for (let i = 0; i < n; i++) for (let d = 0; d < 3; d++) x[3 * i + d] = coords[i][d];
    const g = new Float64Array(N);
    const freeze = (v) => { for (let i = 0; i < n; i++) if (frozen[i]) v[3 * i] = v[3 * i + 1] = v[3 * i + 2] = 0; };
    let E = energyGrad(ff, x, g); freeze(g);
    const M = 8, S = [], Y = [], RHO = [];
    const dir = new Float64Array(N), xn = new Float64Array(N), gn = new Float64Array(N);
    const alpha = new Float64Array(M);
    for (let it = 0; it < maxIter; it++) {
      let gmax = 0; for (let t = 0; t < N; t++) gmax = Math.max(gmax, Math.abs(g[t]));
      if (gmax < 0.02) break;
      // two-loop recursion
      for (let t = 0; t < N; t++) dir[t] = -g[t];
      for (let q = S.length - 1; q >= 0; q--) {
        let a = 0; for (let t = 0; t < N; t++) a += S[q][t] * dir[t];
        a *= RHO[q]; alpha[q] = a;
        for (let t = 0; t < N; t++) dir[t] -= a * Y[q][t];
      }
      if (S.length) {
        const q = S.length - 1; let sy = 0, yy = 0;
        for (let t = 0; t < N; t++) { sy += S[q][t] * Y[q][t]; yy += Y[q][t] * Y[q][t]; }
        const gam = sy / yy; for (let t = 0; t < N; t++) dir[t] *= gam;
      }
      for (let q = 0; q < S.length; q++) {
        let b = 0; for (let t = 0; t < N; t++) b += Y[q][t] * dir[t];
        b *= RHO[q];
        for (let t = 0; t < N; t++) dir[t] += S[q][t] * (alpha[q] - b);
      }
      freeze(dir);
      let slope = 0; for (let t = 0; t < N; t++) slope += dir[t] * g[t];
      if (!(slope < 0)) {                     // not a descent direction: reset to steepest descent
        S.length = Y.length = RHO.length = 0;
        for (let t = 0; t < N; t++) dir[t] = -g[t];
        freeze(dir);
        slope = 0; for (let t = 0; t < N; t++) slope += dir[t] * g[t];
      }
      // cap the per-atom displacement (0.3 A; 0.1 A on the very first steepest step)
      let dmax = 0;
      for (let i = 0; i < n; i++) dmax = Math.max(dmax, Math.hypot(dir[3 * i], dir[3 * i + 1], dir[3 * i + 2]));
      const cap = S.length ? 0.3 : 0.1;
      let step = dmax > cap ? cap / dmax : 1;
      let En = 0, ok = false;
      for (let ls = 0; ls < 30; ls++) {
        for (let t = 0; t < N; t++) xn[t] = x[t] + step * dir[t];
        En = energyGrad(ff, xn, gn); freeze(gn);
        if (En <= E + 1e-4 * step * slope) { ok = true; break; }
        step *= 0.5;
      }
      if (!ok) {
        if (S.length) { S.length = Y.length = RHO.length = 0; continue; }
        break;
      }
      const s = new Float64Array(N), y = new Float64Array(N);
      let sy = 0;
      for (let t = 0; t < N; t++) { s[t] = xn[t] - x[t]; y[t] = gn[t] - g[t]; sy += s[t] * y[t]; }
      if (sy > 1e-12) {
        S.push(s); Y.push(y); RHO.push(1 / sy);
        if (S.length > M) { S.shift(); Y.shift(); RHO.shift(); }
      }
      const dE = E - En;
      x.set(xn); g.set(gn); E = En;
      if (dE < 1e-7) break;
    }
    const out = [];
    for (let i = 0; i < n; i++) out.push([x[3 * i], x[3 * i + 1], x[3 * i + 2]]);
    return out;
  }

  // ================================================================ bonds & rotation
  function degree(mol, adj, i) { return adj[i].length + (mol.atoms[i].hCount || 0); }

  const REASONS = {
    double: "This is a double bond. Twisting it would break the π bond (the side-by-side p orbitals), so it can't rotate.",
    triple: "Triple bond: both ends are linear, so there's nothing to rotate.",
    ring: "This bond is in a ring. Spinning it would tear the ring open, so it can only flex (like a chair flip), not rotate freely.",
    sp: 'One end is linear (sp), so rotating it changes nothing.'
  };

  function selectableBonds(mol) {
    const adj = adjacency(mol), out = [];
    mol.bonds.forEach((b, bi) => {
      if (isH(mol, b.a) || isH(mol, b.b)) return;
      if (degree(mol, adj, b.a) < 2 || degree(mol, adj, b.b) < 2) return;
      const order = b.order || 1;
      let reason = null;
      if (order === 2) reason = REASONS.double;
      else if (order === 3) reason = REASONS.triple;
      else if (b.inRing || bondInCycle(adj, b.a, b.b)) reason = REASONS.ring;
      else if (hybridOf(mol, adj, b.a) === 'sp' || hybridOf(mol, adj, b.b) === 'sp') reason = REASONS.sp;
      out.push({ bond: bi, a: Math.min(b.a, b.b), b: Math.max(b.a, b.b), rotatable: reason === null, reason });
    });
    return out;
  }

  function rotatableBonds(mol) { return selectableBonds(mol).filter(b => b.rotatable).map(b => b.bond); }

  function sideOf(mol, front, back) {
    const adj = adjacency(mol), seen = new Set([back]), q = [back];
    while (q.length) {
      const x = q.shift();
      for (const e of adj[x]) {
        if (e.n === front) { if (x !== back) return []; continue; }   // ring: front reachable another way
        if (!seen.has(e.n)) { seen.add(e.n); q.push(e.n); }
      }
    }
    return Array.from(seen).sort((a, b) => a - b);
  }

  function rotateAbout(mol, coords, front, back, deltaDeg) {
    const out = coords.map(p => p.slice());
    const side = sideOf(mol, front, back);
    if (!side.length || !deltaDeg) return out;
    const o = coords[front], k = unit(sub(coords[back], o)), ang = deltaDeg * DEG;
    side.forEach(i => { out[i] = add(o, rotVec(sub(coords[i], o), k, ang)); });
    return out;
  }

  // IUPAC signed dihedral in (-180, 180]
  function dihedral(coords, i, j, k, l) {
    const b1 = sub(coords[j], coords[i]), b2 = sub(coords[k], coords[j]), b3 = sub(coords[l], coords[k]);
    const n1 = cross(b1, b2), n2 = cross(b2, b3);
    const d = Math.atan2(len(b2) * dot(b1, n2), dot(n1, n2)) / DEG;
    return d <= -180 ? 180 : d;
  }

  function neighborsExcept(mol, adj, i, ex) { return nbrs(adj, i).filter(x => x !== ex); }

  // All substituents on `at` that share the top score (ties: two CH3 on one carbon).
  function candSubs(mol, adj, at, other) {
    const cand = neighborsExcept(mol, adj, at, other);
    let bs = -Infinity;
    const sc = cand.map(x => { const v = groupOf(mol, x, at).score || 0; if (v > bs) bs = v; return v; });
    return cand.filter((x, k) => sc[k] > bs - 1e-9);
  }
  function bestSub(mol, adj, at, other) {
    const c = candSubs(mol, adj, at, other);
    return c.length ? c[0] : -1;
  }

  // Reference pair (first of each tie list) plus the tie lists. Ties are only kept when
  // neither end is a symmetric rotor (for CH3 or C(CH3)3 ends every position is equivalent).
  function priorityPair(mol, front, back) {
    const adj = adjacency(mol);
    const ti = candSubs(mol, adj, front, back), tl = candSubs(mol, adj, back, front);
    const sym = symmetricEnd(mol, front, back) || symmetricEnd(mol, back, front);
    return {
      i: ti.length ? ti[0] : -1, l: tl.length ? tl[0] : -1,
      ti: sym ? ti.slice(0, 1) : ti, tl: sym ? tl.slice(0, 1) : tl
    };
  }

  // Among tied pairs, the one that gives the most specific name: anti, then totally
  // eclipsed, then gauche, then eclipsed; otherwise the reference pair. Equal ranks prefer
  // the positive angle, then the reference order.
  function pairRank(d) {
    const a = Math.abs(d);
    if (a >= 180 - TOL) return 0;
    if (a <= TOL) return 1;
    if (Math.abs(a - 60) <= TOL) return 2;
    if (Math.abs(a - 120) <= TOL) return 3;
    return 4;
  }
  function pickPair(list) {
    let best = null;
    list.forEach((p, k) => {
      const r = pairRank(p.deg);
      const better = !best || r < best.r || (r === best.r && r < 4 && p.deg > 0 && best.p.deg <= 0);
      if (better) best = { p, r, k };
    });
    return best ? best.p : null;
  }

  function dihedralOfPriorityGroups(mol, coords, front, back) {
    const p = priorityPair(mol, front, back);
    if (p.i < 0 || p.l < 0) return { deg: 0, i: p.i, l: p.l };
    const list = [];
    p.ti.forEach(i => p.tl.forEach(l => {
      let d = dihedral(coords, i, front, back, l);
      if (!isFinite(d)) d = 0;
      list.push({ deg: d, i, l });
    }));
    const best = pickPair(list);
    return { deg: best.deg, i: best.i, l: best.l };
  }

  // Rotate so the displayed priority dihedral reads targetDeg. With tied groups several
  // rotations qualify; the first one whose displayed pair then reads the target (or its
  // mirror) wins, else the reference pair is used.
  function setPriorityDihedral(mol, coords, front, back, targetDeg) {
    const p = priorityPair(mol, front, back);
    if (p.i < 0 || p.l < 0) return coords.map(q => q.slice());
    let mirror = null;
    for (const t of [targetDeg, -targetDeg]) for (const i of p.ti) for (const l of p.tl) {
      const cur = dihedral(coords, i, front, back, l);
      const out = rotateAbout(mol, coords, front, back, wrap180(t - cur));
      const shown = dihedralOfPriorityGroups(mol, out, front, back).deg;
      if (Math.abs(wrap180(shown - targetDeg)) < 1) return out;
      if (!mirror && Math.abs(Math.abs(shown) - Math.abs(targetDeg)) < 1) mirror = out;
    }
    if (mirror) return mirror;
    const cur = dihedral(coords, p.i, front, back, p.l);
    return rotateAbout(mol, coords, front, back, wrap180(targetDeg - cur));
  }

  // Rotate so the reference pair (the energy plot's x axis) sits at phiDeg.
  function setReferencePhi(mol, coords, front, back, phiDeg) {
    const p = priorityPair(mol, front, back);
    if (p.i < 0 || p.l < 0) return coords.map(q => q.slice());
    const cur = dihedral(coords, p.i, front, back, p.l);
    return rotateAbout(mol, coords, front, back, wrap180(phiDeg - cur));
  }

  // ================================================================ Newman geometry
  const LP_HALF = 54.7356;   // half the lone-pair / lone-pair opening (tetrahedral)

  function lonePairs(mol, adj, coords, at) {
    const el = mol.atoms[at].el;
    if (el !== 'O' && el !== 'N' && el !== 'S') return [];
    if (hybridOf(mol, adj, at) !== 'sp3') return [];
    const nb = nbrs(adj, at), u = nb.map(x => unit(sub(coords[x], coords[at])));
    if ((el === 'O' || el === 'S') && nb.length === 2) {
      const d = unit(scl(add(u[0], u[1]), -1)), nn = unit(cross(u[0], u[1]));
      const c = Math.cos(LP_HALF * DEG), s = Math.sin(LP_HALF * DEG);
      return [add(scl(d, c), scl(nn, s)), sub(scl(d, c), scl(nn, s))];
    }
    if (el === 'N' && nb.length === 3) return [unit(scl(add(add(u[0], u[1]), u[2]), -1))];
    return [];
  }

  function newmanAngles(mol, coords, front, back, upVec) {
    const adj = adjacency(mol);
    const w = unit(sub(coords[back], coords[front]));
    let u = null, upDeg = 0;
    if (upVec) {
      u = unit(perp(upVec, w));
      if (len(u) < 0.5) u = null;
    }
    if (!u) {
      const pi = bestSub(mol, adj, front, back);
      if (pi >= 0) u = unit(perp(sub(coords[pi], coords[front]), w));
      if (!u || len(u) < 0.5) u = anyPerp(w);
    }
    const rr = cross(w, u);
    const ang = (v) => wrap360(Math.atan2(dot(v, rr), dot(v, u)) / DEG);
    const frontSubs = neighborsExcept(mol, adj, front, back).map(x => ({ atom: x, angle: ang(perp(sub(coords[x], coords[front]), w)) }));
    const backSubs = neighborsExcept(mol, adj, back, front).map(x => ({ atom: x, angle: ang(perp(sub(coords[x], coords[back]), w)) }));
    lonePairs(mol, adj, coords, front).forEach(v => frontSubs.push({ atom: null, lp: true, angle: ang(perp(v, w)) }));
    lonePairs(mol, adj, coords, back).forEach(v => backSubs.push({ atom: null, lp: true, angle: ang(perp(v, w)) }));
    if (upVec) {
      const pi = bestSub(mol, adj, front, back);
      const f = frontSubs.find(s => s.atom === pi);
      upDeg = f ? f.angle : 0;
    }
    return { frontSubs, backSubs, upDeg };
  }

  // ================================================================ energy model (SPEC 4.8)
  const ALKYL = ['H', 'Me', 'iPr', 'tBu'];
  const ECL_ALK = {
    'H|H': 4.0, 'H|Me': 6.0, 'H|iPr': 7.0, 'H|tBu': 8.0,
    'Me|Me': 11.0, 'Me|iPr': 13.0, 'Me|tBu': 16.0,
    'iPr|iPr': 16.0, 'iPr|tBu': 20.0, 'tBu|tBu': 28.0
  };
  const ECL_H_X = { F: 4.0, Cl: 5.0, Br: 5.0, I: 5.0, O: 4.5, N: 5.0, S: 5.0 };
  const ECL_ME_X = { F: 8, Cl: 12, Br: 13, I: 14, O: 10, N: 10, S: 12 };
  const ECL_HOMO = { F: 10, Cl: 16, Br: 19, I: 22, O: 12, N: 12, S: 16 };
  const G_ALK = { 'Me|Me': 3.8, 'Me|iPr': 5.0, 'Me|tBu': 11, 'iPr|iPr': 7.0, 'iPr|tBu': 14, 'tBu|tBu': 25 };
  const G_ME_X = { F: 0.8, Cl: 2.0, Br: 2.5, I: 3.0, O: 1.5, N: 2.0, S: 2.5 };
  const ALK_PLUS_ECL = { Me: 0, iPr: 2, tBu: 4 };
  const ALK_PLUS_G = { Me: 0, iPr: 1, tBu: 3 };
  const alkKey = (a, b) => (ALKYL.indexOf(a) <= ALKYL.indexOf(b) ? a + '|' + b : b + '|' + a);

  function pairEcl(c1, c2) {
    const a1 = ALKYL.includes(c1), a2 = ALKYL.includes(c2);
    if (a1 && a2) return ECL_ALK[alkKey(c1, c2)];
    if (a1 || a2) {
      const alk = a1 ? c1 : c2, x = a1 ? c2 : c1;
      if (alk === 'H') return ECL_H_X[x];
      return ECL_ME_X[x] + ALK_PLUS_ECL[alk];
    }
    return (ECL_HOMO[c1] + ECL_HOMO[c2]) / 2;
  }
  function pairG(c1, c2, h1, h2) {
    if (c1 === 'H' || c2 === 'H') return 0;
    const a1 = ALKYL.includes(c1), a2 = ALKYL.includes(c2);
    if (a1 && a2) return G_ALK[alkKey(c1, c2)];
    if (a1 || a2) {
      const alk = a1 ? c1 : c2, x = a1 ? c2 : c1;
      return G_ME_X[x] + ALK_PLUS_G[alk];
    }
    if (c1 === c2) {
      switch (c1) {
        case 'F': return -2.4;
        case 'Cl': return 4.6;
        case 'Br': return 6.3;
        case 'I': return 8.0;
        case 'S': return 3.0;
        case 'N': return 1.0;
        case 'O': return (h1 || h2) ? -3.0 : 1.0;
      }
    }
    const pair = [c1, c2].sort().join('|');
    const oH = (c1 === 'O' && h1) || (c2 === 'O' && h2);
    if (pair === 'F|O') return oH ? -2.0 : -1.0;
    if (pair === 'N|O') return (h1 || h2) ? -2.0 : 1.0;
    return 1.5;
  }

  const DISPLAY = { H: 'H', Me: 'CH3', iPr: 'CHR2', tBu: 'CR3' };
  function clsDisplay(c) { return DISPLAY[c] || c; }

  function ecl(theta) { const a = Math.abs(theta); return a <= 60 ? Math.pow(Math.cos(1.5 * a * DEG), 2) : 0; }
  function gau(theta) { const a = Math.abs(theta); return a <= 120 ? Math.pow(Math.sin(1.5 * a * DEG), 2) : 0; }

  // Everything the energy model needs for one bond: substituent classes + Newman angles.
  // Two kinds: 'alkyl' (sp3-sp3, additive pair model) and 'allylic' (sp3 next to C=C or C=O).
  function energyContext(mol, coords, front, back) {
    const b = mol.bonds[bondIndex(mol, front, back)];
    if (!b || (b.order || 1) !== 1) return null;
    const adj = adjacency(mol);
    if (b.inRing || bondInCycle(adj, front, back)) return null;
    const hf = hybridOf(mol, adj, front), hb = hybridOf(mol, adj, back);
    if (hf === 'sp' || hb === 'sp') return null;
    const na = newmanAngles(mol, coords, front, back);
    const fr = na.frontSubs.filter(s => !s.lp).map(s => ({ atom: s.atom, angle: s.angle, g: groupOf(mol, s.atom, front) }));
    const bk = na.backSubs.filter(s => !s.lp).map(s => ({ atom: s.atom, angle: s.angle, g: groupOf(mol, s.atom, back) }));
    const pr = priorityPair(mol, front, back);
    const fi = fr.find(s => s.atom === pr.i), bl = bk.find(s => s.atom === pr.l);
    if (!fi || !bl) return null;
    const cur = wrap180(bl.angle - fi.angle);
    const angOf = (x) => { const s = fr.find(q => q.atom === x) || bk.find(q => q.atom === x); return s ? s.angle : 0; };
    // tied priority pairs as offsets from the reference pair (for naming)
    const ties = [];
    pr.ti.forEach(i => pr.tl.forEach(l => ties.push({ i, l, off: wrap180(angOf(l) - angOf(i) - cur) })));

    if (hf === 'sp2' || hb === 'sp2') {
      if (hf === hb) return null;                     // conjugated sp2-sp2 single bond: not modelled
      const s2 = hf === 'sp2' ? front : back, s3 = s2 === front ? back : front;
      const e2 = adj[s2].find(e => e.n !== s3 && e.order === 2);
      if (!e2) return null;
      const Y = e2.n, yEl = mol.atoms[Y].el, Ya = angOf(Y);
      const xs = (s3 === front ? fr : bk);
      const cc = yEl === 'C';
      const terms = xs.map(x => {
        const heavy = x.g.cls !== 'H';
        const d = x.angle - Ya;
        return {
          atom: x.atom, cls: x.g.cls, heavy,
          rel: (s3 === back ? d : -d) - cur,
          c: cc ? (heavy ? 0.567 : 1.367) : (heavy ? 0.35 : 0.8),
          // heavy X: syn (eclipsing C=C) +0.8 above skew; for C=O syn is the minimum, skew +3.8
          q: heavy ? (cc ? 0.533 : -2.533) : 0, r: heavy ? (cc ? 0.267 : 2.533) : 0
        };
      });
      const sym = xs.length >= 2 && xs.every(x => x.g.cls === xs[0].g.cls);
      const xPri = terms.find(t => t.atom === (s3 === front ? pr.i : pr.l)) || terms[0];
      return { kind: 'allylic', terms, cur, pr, fr, bk, Y, yEl, s2, s3, sym, xPri, ties };
    }
    const pairs = [];
    fr.forEach(f => bk.forEach(k => {
      pairs.push({
        rel: k.angle - f.angle - cur,       // pair dihedral = phi + rel
        E: pairEcl(f.g.cls, k.g.cls), G: pairG(f.g.cls, k.g.cls, f.g.hasH, k.g.hasH),
        f: f.g, k: k.g
      });
    }));
    return { kind: 'alkyl', pairs, cur, pr, fr, bk, ties };
  }
  function rawEnergy(ctx, phi) {
    let E = 0;
    if (ctx.kind === 'allylic') {
      for (const t of ctx.terms) {
        const th = (phi + t.rel) * DEG;
        E += t.c * (1 - Math.cos(3 * th)) + (t.heavy ? t.q * Math.cos(2 * th) + t.r : 0);
      }
      return E;
    }
    for (const p of ctx.pairs) {
      const th = wrap180(phi + p.rel);
      E += p.E * ecl(th) + p.G * gau(th);
    }
    return E;
  }
  function gridMin(ctx) {
    let mn = Infinity;
    for (let d = 0; d <= 360; d += 2) mn = Math.min(mn, rawEnergy(ctx, d));
    return mn;
  }

  // Name at a reference-pair angle phi (the energy plot's x axis).
  function nameAt(mol, ctx, front, back, phi) {
    if (ctx.kind === 'allylic') return allylicName(ctx, phi);
    const list = ctx.ties.map(t => ({ deg: wrap180(phi + t.off) }));
    const best = pickPair(list) || { deg: phi };
    return angleName(mol, front, back, best.deg);
  }

  const ALLYL_TEXT = {
    syn: ['Syn (eclipses the C=Y)', 'syn'], 'allyl-skew': ['Skew', 'skew'], eclipsed: ['Eclipsed', 'eclipsed'],
    between: ['In between', 'in between'], bisected: ['Bisected', 'bisected']
  };
  function allylicName(ctx, phi) {
    const y = ctx.yEl;
    if (ctx.sym) {
      const a = Math.abs(wrap180(phi + ctx.terms[0].rel)), r = a % 120;
      if (near(r, 0) || near(r, 120)) return { key: 'db-eclipsed', text: 'C=' + y + ' eclipsed', name: 'C=' + y + ' eclipsed' };
      if (near(r, 60)) return { key: 'bisected', text: 'Bisected', name: 'bisected' };
      return { key: 'skew', text: TEXT.skew, name: SHORT.skew };
    }
    const a = Math.abs(wrap180(phi + ctx.xPri.rel));
    let key;
    if (near(a, 0)) key = 'syn';
    else if (a >= 100 && a <= 150) key = 'allyl-skew';
    else if (near(a, 60) || near(a, 180)) key = 'eclipsed';
    else key = 'between';
    const t = ALLYL_TEXT[key];
    return { key, text: t[0].replace('C=Y', 'C=' + y), name: t[1] };
  }

  function bondIndex(mol, a, b) {
    for (let i = 0; i < mol.bonds.length; i++) {
      const q = mol.bonds[i];
      if ((q.a === a && q.b === b) || (q.a === b && q.b === a)) return i;
    }
    return -1;
  }

  function torsionProfile(mol, coords, front, back) {
    const ctx = energyContext(mol, coords, front, back);
    if (!ctx) return null;
    const raw = [];
    for (let d = 0; d <= 360; d += 2) raw.push(rawEnergy(ctx, d));
    const mn = Math.min.apply(null, raw);
    const points = raw.map((e, k) => ({ deg: 2 * k, kJ: e - mn }));
    const curDeg = wrap360(ctx.cur);
    const current = { deg: curDeg, kJ: Math.max(0, rawEnergy(ctx, curDeg) - mn) };
    let max = 0; points.forEach(p => { if (p.kJ > max) max = p.kJ; });

    let terms = [], note = null;
    if (ctx.kind === 'allylic') {
      const yl = 'C=' + ctx.yEl;
      note = 'Allylic teaching model (approximate kJ/mol): a C–H or C–C bond on the sp3 carbon prefers to eclipse the ' + yl +
        ' bond, the opposite of the ethane rule. ' + (ctx.yEl === 'C'
        ? 'Propene barrier ≈ 8.2; for a CH3 group, syn ≈ 0.8 above skew, maxima ≈ 7.'
        : 'Barrier ≈ 4.8 for H; a CH3 group prefers to eclipse the ' + yl + ' (skew ≈ 3.8 higher).');
    } else {
      // distinct class pairs actually used
      const seen = new Map();
      ctx.pairs.forEach(p => {
        const c1 = p.f.cls, c2 = p.k.cls;
        const oh = (c1 === 'O' || c2 === 'O' || c1 === 'N' || c2 === 'N') && ALKYL.indexOf(c1) < 0 && ALKYL.indexOf(c2) < 0 && !!(p.f.hasH || p.k.hasH);
        const key = [c1, c2].sort().join('|') + (oh ? '*' : '');
        if (!seen.has(key)) seen.set(key, { front: c1, back: c2, ecl: p.E, gauche: p.G, label: clsDisplay(c1) + '/' + clsDisplay(c2), oh });
      });
      terms = Array.from(seen.values());
    }

    // stationary points on the circular 2-degree grid (window +-20 deg to ignore wiggles)
    const U = raw.length - 1, W = 10, stationary = [];
    const at = (k) => points[((k % U) + U) % U].kJ;
    for (let k = 0; k < U; k++) {
      const e = at(k);
      let isMax = true, isMin = true;
      for (let o = 1; o <= W; o++) {
        const l = at(k - o), r = at(k + o);
        if (!(e >= l && e >= r)) isMax = false;
        if (!(e <= l && e <= r)) isMin = false;
      }
      if (isMax && isMin) continue;   // flat
      if (!isMax && !isMin) continue;
      if (isMax && e - Math.min(at(k - W), at(k + W)) < 0.05) continue;
      if (isMin && Math.max(at(k - W), at(k + W)) - e < 0.05) continue;
      const kind = isMax ? 'max' : 'min';
      if (stationary.some(s => s.kind === kind && Math.abs(wrap180(s.deg - 2 * k)) < 2 * W)) continue;
      stationary.push({ deg: 2 * k, kJ: e, kind, name: nameAt(mol, ctx, front, back, 2 * k).name });
    }
    return {
      points, current, priority: { i: ctx.pr.i, l: ctx.pr.l }, kind: ctx.kind, note,
      terms, max, barrier: max, stationary
    };
  }

  // Energy (relative to the profile minimum) at a reference-pair angle, or at the current
  // geometry when priorityDeg is null/undefined.
  function energyAt(mol, coords, front, back, priorityDeg) {
    const ctx = energyContext(mol, coords, front, back);
    if (!ctx) return 0;
    const phi = priorityDeg == null ? ctx.cur : priorityDeg;
    return Math.max(0, rawEnergy(ctx, phi) - gridMin(ctx));
  }

  // Current reference-pair angle (0-360) and energy, for the plot marker.
  function energyState(mol, coords, front, back) {
    const ctx = energyContext(mol, coords, front, back);
    if (!ctx) return null;
    return { deg: wrap360(ctx.cur), kJ: Math.max(0, rawEnergy(ctx, ctx.cur) - gridMin(ctx)) };
  }

  // ================================================================ naming (SPEC 4.9)
  const TOL = 15;
  const near = (a, t) => Math.abs(a - t) <= TOL;

  // Symmetric rotor: all substituents on this end share one energy class (deviation 5)
  function symmetricEnd(mol, at, other) {
    const adj = adjacency(mol);
    const cls = neighborsExcept(mol, adj, at, other).map(x => groupOf(mol, x, at).cls);
    // needs >= 2 substituents: a lone group (the H of an OH, the CH3 of CH3-S-S) is not a rotor
    return cls.length >= 2 && cls.every(c => c === cls[0]);
  }

  const TEXT = {
    anti: 'Anti (staggered)', gauche: 'Gauche (staggered)', eclipsed: 'Eclipsed',
    'totally-eclipsed': 'Totally eclipsed', staggered: 'Staggered', skew: 'Skew (in between)'
  };
  const SHORT = {
    anti: 'anti', gauche: 'gauche', eclipsed: 'eclipsed', 'totally-eclipsed': 'totally eclipsed',
    staggered: 'staggered', skew: 'skew'
  };

  // Name for an open-chain single bond at a given priority dihedral (deg, any range)
  function angleName(mol, front, back, deg) {
    const a = Math.abs(wrap180(deg));
    const pr = priorityPair(mol, front, back);
    const ethaneLike = symmetricEnd(mol, front, back) || symmetricEnd(mol, back, front);
    let key;
    if (ethaneLike) {
      const r = a % 120;   // staggered every 120 deg offset by 60, eclipsed at multiples of 120
      if (near(r, 0) || near(r, 120)) key = 'eclipsed';
      else if (near(r, 60)) key = 'staggered';
      else key = 'skew';
    } else {
      const bothHeavy = pr.i >= 0 && pr.l >= 0 && !isH(mol, pr.i) && !isH(mol, pr.l);
      if (near(a, 180)) key = 'anti';
      else if (near(a, 60)) key = 'gauche';
      else if (near(a, 0)) key = bothHeavy ? 'totally-eclipsed' : 'eclipsed';
      else if (near(a, 120)) key = 'eclipsed';
      else key = 'skew';
    }
    return { key, text: TEXT[key], name: SHORT[key] };
  }

  function ringConformation(mol, coords, ring) {
    const k = ring.length;
    const tau = [];
    for (let j = 0; j < k; j++) tau.push(dihedral(coords, ring[j], ring[(j + 1) % k], ring[(j + 2) % k], ring[(j + 3) % k]));
    const abs = tau.map(Math.abs), mx = Math.max.apply(null, abs);
    if (mx < 10) return 'planar';
    if (k === 6) {
      const alt = tau.every((t, j) => Math.sign(t) !== Math.sign(tau[(j + 1) % k]));
      if (alt && abs.every(a => a >= 35 && a <= 75)) return 'chair';
      const small = []; abs.forEach((a, j) => { if (a < 15) small.push(j); });
      if (small.length >= 2) {
        const opposite = small.some(p => small.includes((p + 3) % 6));
        return opposite ? 'boat' : 'half-chair';
      }
      if (small.length === 1) return 'half-chair';
      if (abs.every(a => a >= 15 && a <= 75)) return 'twist-boat';
      return 'puckered';
    }
    if (k === 5) return 'envelope';
    return 'puckered';
  }

  function conformationName(mol, coords, front, back) {
    const bi = bondIndex(mol, front, back), b = mol.bonds[bi];
    const adj = adjacency(mol);
    const pd = () => { const d = dihedralOfPriorityGroups(mol, coords, front, back).deg; return isFinite(d) ? d : 0; };
    if (!b) return { key: 'fixed', text: 'Not a bond', dihedral: 0 };
    const order = b.order || 1;
    if (order === 2) {
      const pp = priorityPair(mol, front, back);
      if (pp.i >= 0 && pp.l >= 0 && !isH(mol, pp.i) && !isH(mol, pp.l)) {
        const d2 = Math.abs(pd());
        if (d2 > 150) return { key: 'fixed', text: 'Trans (double bond)', dihedral: pd() };
        if (d2 < 30) return { key: 'fixed', text: 'Cis (double bond)', dihedral: pd() };
      }
      return { key: 'fixed', text: 'Double bond (planar, fixed)', dihedral: pd() };
    }
    if (order === 3) return { key: 'fixed', text: 'Triple bond (linear, fixed)', dihedral: 0 };
    if (b.inRing || bondInCycle(adj, front, back)) {
      const rings = ringsOf(mol, adj).filter(r => r.includes(front) && r.includes(back));
      rings.sort((p, q) => p.length - q.length);
      const ring = rings[0];
      if (!ring) return { key: 'ring', text: 'Ring bond', dihedral: 0, ring: true };
      const k = ring.length;
      const conf = ringConformation(mol, coords, ring);
      if (k <= 4) {
        // 3- and 4-rings: the ring torsion is ~0 by construction; report the cis H-C-C-H
        // (or X-C-C-Y) dihedral, which shows the eclipsing strain.
        const fs = nbrs(adj, front).filter(x => !ring.includes(x)), bs = nbrs(adj, back).filter(x => !ring.includes(x));
        let best = null;
        fs.forEach(x => bs.forEach(y => { const d = dihedral(coords, x, front, back, y); if (!best || Math.abs(d) < Math.abs(best.d)) best = { d, x, y }; }));
        const tor = best ? best.d : 0;
        const key = Math.abs(tor) <= TOL ? 'ring-eclipsed' : 'ring';
        return {
          key, ring: true, size: k, dihedral: tor, torsionAtoms: best ? [best.x, front, back, best.y] : null,
          text: (k === 3 ? 'Eclipsed (flat 3-ring)' : Math.abs(tor) <= 30 ? 'Nearly eclipsed (4-ring)' : 'Ring bond (4-ring)')
        };
      }
      // ring torsion through this bond: prev(front) - front - back - next(back) along the ring
      const fi = ring.indexOf(front), bi2 = ring.indexOf(back);
      const dirn = ((fi + 1) % k === bi2) ? 1 : -1;
      const pf = ring[(fi - dirn + k) % k], nb = ring[(bi2 + dirn + k) % k];
      const tor = dihedral(coords, pf, front, back, nb);
      const out = { dihedral: tor, ring: true, size: k, torsionAtoms: [pf, front, back, nb] };
      if (conf === 'chair' || conf === 'boat' || conf === 'twist-boat') {
        return Object.assign(out, { key: conf, text: conf.charAt(0).toUpperCase() + conf.slice(1) });
      }
      return Object.assign(out, { key: 'ring', text: 'Ring bond (' + conf + ')' });
    }
    const hf = hybridOf(mol, adj, front), hb = hybridOf(mol, adj, back);
    if (hf === 'sp' || hb === 'sp') {
      return { key: 'fixed', text: 'Linear end (fixed)', dihedral: pd() };
    }
    const d = pd();
    if (hf === 'sp2' || hb === 'sp2') {
      if (hf === hb) {
        // conjugated single bond (butadiene C2-C3): s-trans / s-cis
        const y1 = adj[front].find(e => e.n !== back && e.order === 2), y2 = adj[back].find(e => e.n !== front && e.order === 2);
        if (y1 && y2) {
          const t = Math.abs(dihedral(coords, y1.n, front, back, y2.n));
          const key = near(t, 180) ? 's-trans' : near(t, 0) ? 's-cis' : 'twisted';
          return { key, text: key === 's-trans' ? 's-trans' : key === 's-cis' ? 's-cis' : 'Twisted', dihedral: d, conj: true };
        }
        return { key: 'skew', text: TEXT.skew, dihedral: d, conj: true };
      }
      const ctx = energyContext(mol, coords, front, back);
      if (ctx) { const nm = allylicName(ctx, ctx.cur); return { key: nm.key, text: nm.text, dihedral: d, allylic: true }; }
    }
    const nm = angleName(mol, front, back, d);
    return { key: nm.key, text: nm.text, dihedral: d };
  }

  // Axial / equatorial for every non-ring substituent of a ring atom (sp3 rings only).
  // Returns a Map atom -> 'ax' | 'eq'.
  function axialEquatorial(mol, coords) {
    const adj = adjacency(mol), out = new Map();
    ringsOf(mol, adj).forEach(ring => {
      if (ring.length < 5) return;
      const cen = ring.reduce((s, a) => add(s, coords[a]), [0, 0, 0]).map(v => v / ring.length);
      let nn = [0, 0, 0];
      for (let j = 0; j < ring.length; j++) nn = add(nn, cross(sub(coords[ring[j]], cen), sub(coords[ring[(j + 1) % ring.length]], cen)));
      nn = unit(nn);
      ring.forEach(a => {
        if (hybridOf(mol, adj, a) !== 'sp3') return;
        nbrs(adj, a).filter(x => !ring.includes(x)).forEach(x => {
          const c = Math.abs(dot(unit(sub(coords[x], coords[a])), nn));
          out.set(x, c > 0.7 ? 'ax' : 'eq');
        });
      });
    });
    return out;
  }

  // ================================================================ exports
  const PAIR_TABLE = {
    ecl: {
      alkyl: ECL_ALK, 'H-X': ECL_H_X, 'Me-X': ECL_ME_X, 'X-X homo (mixed = mean)': ECL_HOMO,
      note: 'iPr-X = Me-X + 2, tBu-X = Me-X + 4'
    },
    gauche: {
      alkyl: G_ALK, 'Me-X': G_ME_X,
      hetero: { 'F|F': -2.4, 'Cl|Cl': 4.6, 'Br|Br': 6.3, 'I|I': 8.0, 'S|S': 3.0, 'N|N': 1.0, 'O|O': 1.0, 'O|O (OH)': -3.0, 'F|O': -1.0, 'F|O (OH)': -2.0, 'N|O': 1.0, 'N|O (with H)': -2.0, other: 1.5 },
      note: 'Any pair with H = 0. iPr-X = Me-X + 1, tBu-X = Me-X + 3.'
    },
    note: 'Approximate teaching values (kJ/mol), additive pair model.'
  };

  const api = {
    build3D, relax, selectableBonds, rotatableBonds, sideOf, rotateAbout, dihedral,
    priorityPair, dihedralOfPriorityGroups, setPriorityDihedral, setReferencePhi, newmanAngles,
    torsionProfile, energyAt, energyState, conformationName, ringConformation, axialEquatorial, PAIR_TABLE,
    // internals exposed for tests / debugging (not part of the page contract)
    _internal: { setupFF, energyGrad, embed, mulberry32, pairEcl, pairG, bondR0, SEED }
  };
  root.NNGeom = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

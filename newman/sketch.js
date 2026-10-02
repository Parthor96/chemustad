/*
 * newman/sketch.js  —  window.NNSketch (Newman Navigator, skeletal-structure sketcher)
 *
 * The student draws a line (skeletal) structure with a mouse, a finger, or a pen.
 * Freehand strokes are cleaned up into straight bonds (resample → Ramer–Douglas–Peucker →
 * merge near-collinear turns → drop tiny segments → snap ends onto existing atoms).
 * A click-to-place "chain" mode builds zigzags one atom at a time.
 * Contract: newman/SPEC.md section 5. Pointer Events only (mouse, touch, Apple Pencil).
 *
 * The graph lives in CSS pixels: vertices [{x, y, el, charge}], edges [{a, b, order}].
 * getMolecule() hands it to NNChem.fromGraph(..., 'sketch'), which adds implicit H by valence.
 *
 * Resolved details (signatures unchanged):
 *  - Bond-order cycling skips an order that would break a valence (e.g. 2 → 3 on a carbon
 *    that already has three other bonds goes back to 1); the bond flashes red for 300 ms
 *    so the student sees the triple bond was refused. Without the skip a bond could get
 *    stuck at double with no way back except Undo.
 *  - In chain mode a tap on the active vertex itself opens the element picker (a long
 *    press anywhere on a vertex does too, in both modes).
 *  - In chain mode a new point that lands on an existing atom bonds to that atom
 *    (ring closure) instead of stacking a second atom on top of it.
 *  - The largest allowed valence is the first entry of NNChem.EL[el].valence (the
 *    "normal" valence), so S is treated as divalent in drawings.
 */
(function (root) {
  'use strict';

  const COL = {
    bg: '#17131d', grid: '#2d2436', bond: '#fbf3e2', dot: '#9d9386', err: '#e53935',
    stroke: 'rgba(255,179,0,.7)', teal: '#14b8a6', mari: '#ffb300', dim: '#9d9386'
  };
  const FONT = 'Inter, -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif';
  const ELEMENTS = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I'];
  const FALLBACK_EL = {
    C: { color2d: '#ffb300', valence: [4] }, N: { color2d: '#6f9bff', valence: [3] },
    O: { color2d: '#ff6b61', valence: [2] }, S: { color2d: '#f2c200', valence: [2] },
    F: { color2d: '#9be15d', valence: [1] }, Cl: { color2d: '#2ec27e', valence: [1] },
    Br: { color2d: '#e0875f', valence: [1] }, I: { color2d: '#b784f0', valence: [1] }
  };
  const EL_NAME = { C: 'Carbon', N: 'Nitrogen', O: 'Oxygen', S: 'Sulfur', F: 'Fluorine', Cl: 'Chlorine', Br: 'Bromine', I: 'Iodine' };

  const chem = () => root.NNChem || null;
  const elInfo = (el) => { const c = chem(); return (c && c.EL && c.EL[el]) || FALLBACK_EL[el] || FALLBACK_EL.C; };
  const maxValence = (el) => elInfo(el).valence[0];

  /* ---------------- geometry helpers ---------------- */
  const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
  function segDist(p, a, b) {
    const dx = b.x - a.x, dy = b.y - a.y, L2 = dx * dx + dy * dy;
    let t = L2 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2 : 0;
    t = Math.max(0, Math.min(1, t));
    return Math.hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
  }
  // Evenly spaced samples along the stroke (spacing in px).
  function resample(pts, step) {
    if (pts.length < 2) return pts.slice();
    const out = [{ x: pts[0].x, y: pts[0].y }];
    let carry = 0;
    for (let i = 1; i < pts.length; i++) {
      let a = pts[i - 1]; const b = pts[i];
      let seg = dist(a, b);
      while (carry + seg >= step) {
        const t = (step - carry) / seg;
        const p = { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t };
        out.push(p); a = p; seg = dist(a, b); carry = 0;
      }
      carry += seg;
    }
    const last = pts[pts.length - 1];
    if (dist(out[out.length - 1], last) > 0.5) out.push({ x: last.x, y: last.y });
    return out;
  }
  // Ramer–Douglas–Peucker simplification.
  function rdp(pts, eps) {
    if (pts.length < 3) return pts.slice();
    let idx = -1, dmax = 0;
    const a = pts[0], b = pts[pts.length - 1];
    for (let i = 1; i < pts.length - 1; i++) {
      const d = segDist(pts[i], a, b);
      if (d > dmax) { dmax = d; idx = i; }
    }
    if (dmax <= eps) return [a, b];
    const left = rdp(pts.slice(0, idx + 1), eps), right = rdp(pts.slice(idx), eps);
    return left.slice(0, -1).concat(right);
  }
  function turnAngle(a, b, c) {
    const a1 = Math.atan2(b.y - a.y, b.x - a.x), a2 = Math.atan2(c.y - b.y, c.x - b.x);
    let d = Math.abs(a2 - a1) * 180 / Math.PI;
    if (d > 180) d = 360 - d;
    return d;
  }
  // Freehand stroke → the corner points of the cleaned-up polyline.
  function cleanStroke(raw) {
    let p = resample(raw, 3);
    p = rdp(p, 7);
    // merge near-collinear interior points
    let changed = true;
    while (changed && p.length > 2) {
      changed = false;
      for (let i = 1; i < p.length - 1; i++) {
        if (turnAngle(p[i - 1], p[i], p[i + 1]) < 20) { p.splice(i, 1); changed = true; break; }
      }
    }
    // drop segments shorter than 14 px by merging their endpoints
    changed = true;
    while (changed && p.length > 1) {
      changed = false;
      for (let i = 1; i < p.length; i++) {
        if (dist(p[i - 1], p[i]) < 14) {
          if (p.length === 2) return [];
          if (i === p.length - 1) p.splice(i - 1, 1); else p.splice(i, 1);
          changed = true; break;
        }
      }
    }
    return p;
  }

  /* ---------------- wedge / dash → chirality ---------------- */
  // V: [{x, y, el}], E: [{a, b, order, st?: 'wedge'|'dash', from?: vertex at the narrow end}].
  // For each atom at the narrow end of a wedge or dash, build 3D directions to its neighbours
  // (screen x right, y up, wedge +z toward the viewer, dash -z), add the implicit H opposite the
  // other three, and return the SMILES-style tag: chiral '@' when det(p1-p0, p2-p0, p3-p0) < 0,
  // the same test newman/geom.js uses when it places the atoms. chiralNbrs uses -1 for the implicit H.
  // Returns { tags: {vertex: {chiral, chiralNbrs}}, skipped: [vertex] } (skipped = not a stereocenter shape).
  function stereoFromDrawing(V, E, maxVal) {
    const tags = {}, skipped = [];
    const centers = new Set();
    E.forEach((e) => { if (e.st && e.order === 1 && (e.from === e.a || e.from === e.b)) centers.add(e.from); });
    centers.forEach((i) => {
      const inc = E.filter((e) => e.a === i || e.b === i);
      if (inc.some((e) => e.order !== 1)) { skipped.push(i); return; }
      const nb = inc.map((e) => (e.a === i ? e.b : e.a));
      const h = Math.max(0, maxVal(V[i].el) - nb.length);
      if (nb.length + h !== 4 || h > 1 || nb.length < 3) { skipped.push(i); return; }
      const vec = inc.map((e, k) => {
        const j = nb[k];
        let dx = V[j].x - V[i].x, dy = -(V[j].y - V[i].y);
        const L = Math.hypot(dx, dy) || 1; dx /= L; dy /= L;
        let z = 0;
        if (e.st && e.from === i) z = e.st === 'wedge' ? 0.9 : -0.9;
        const n = Math.hypot(dx, dy, z);
        return [dx / n, dy / n, z / n];
      });
      if (h === 1) {
        let hv = vec.reduce((s, v) => [s[0] - v[0], s[1] - v[1], s[2] - v[2]], [0, 0, 0]);
        if (Math.hypot(hv[0], hv[1], hv[2]) < 1e-3) { skipped.push(i); return; }
        vec.push(hv);
      }
      const [p0, p1, p2, p3] = vec;
      const a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
      const b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
      const c = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
      const d = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
      if (Math.abs(d) < 1e-4) { skipped.push(i); return; }
      tags[i] = { chiral: d < 0 ? '@' : '@@', chiralNbrs: h === 1 ? nb.concat([-1]) : nb.slice() };
    });
    return { tags, skipped };
  }

  /* ---------------- attach ---------------- */
  function attach(canvas, opts) {
    opts = opts || {};
    const ctx = canvas.getContext('2d');
    const picker = opts.picker || null;
    const bondLength = opts.bondLength || 44;
    let mode = opts.mode === 'chain' ? 'chain' : 'free';
    let bondTool = 'plain'; // 'plain' | 'wedge' | 'dash'

    let V = [];            // vertices {x, y, el, charge}
    let E = [];            // edges {a, b, order}
    const history = [];    // snapshots for undo (max 50)
    let active = null;     // chain mode: index of the vertex new atoms attach to
    let pickerFor = -1;    // vertex index the picker is open for
    let stroke = null;     // raw points while drawing
    let shake = null;      // {edge, until}
    let hover = -1;        // hovered vertex (mouse)
    const listeners = [];
    if (opts.onChange) listeners.push(opts.onChange);

    let W = 1, H = 1, dpr = 1;
    canvas.style.touchAction = 'none';

    /* ----- sizing ----- */
    function resize() {
      dpr = Math.min(root.devicePixelRatio || 1, 2);
      const r = canvas.getBoundingClientRect();
      W = Math.max(1, r.width); H = Math.max(1, r.height);
      canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
      redraw();
    }
    let ro = null;
    if (root.ResizeObserver) { ro = new ResizeObserver(resize); ro.observe(canvas); }

    /* ----- graph helpers ----- */
    const snapshot = () => ({ V: V.map((v) => Object.assign({}, v)), E: E.map((e) => Object.assign({}, e)), active });
    function commitStart() { history.push(snapshot()); if (history.length > 50) history.shift(); }
    function edgeIndex(a, b) { return E.findIndex((e) => (e.a === a && e.b === b) || (e.a === b && e.b === a)); }
    function degree(i) { let s = 0; E.forEach((e) => { if (e.a === i || e.b === i) s += e.order; }); return s; }
    function nbrs(i) { const out = []; E.forEach((e) => { if (e.a === i) out.push(e.b); else if (e.b === i) out.push(e.a); }); return out; }
    function nearestVertex(p, R) {
      let best = -1, bd = R;
      V.forEach((v, i) => { const d = dist(v, p); if (d <= bd) { bd = d; best = i; } });
      return best;
    }
    function nearestEdge(p, R) {
      let best = -1, bd = R;
      E.forEach((e, i) => { const d = segDist(p, V[e.a], V[e.b]); if (d <= bd) { bd = d; best = i; } });
      return best;
    }
    function addVertex(p) { V.push({ x: p.x, y: p.y, el: 'C', charge: 0 }); return V.length - 1; }
    function addEdge(a, b) {
      if (a === b || edgeIndex(a, b) >= 0) return;
      const e = { a, b, order: 1 };
      if (bondTool !== 'plain') { e.st = bondTool; e.from = a; }
      E.push(e);
    }
    function deleteVertex(i) {
      const fix = (x) => (x > i ? x - 1 : x);
      E = E.filter((e) => e.a !== i && e.b !== i).map((e) => {
        const o = { a: fix(e.a), b: fix(e.b), order: e.order };
        if (e.st) { o.st = e.st; o.from = fix(e.from); }
        return o;
      });
      V.splice(i, 1);
      if (active === i) active = null; else if (active != null && active > i) active--;
    }
    const snapR = (type) => (type === 'touch' ? 22 : 16);
    const bondR = (type) => (type === 'touch' ? 16 : 10);

    /* ----- molecule out ----- */
    function pieces() {
      if (!V.length) return 0;
      const seen = new Array(V.length).fill(false);
      let n = 0;
      for (let s = 0; s < V.length; s++) {
        if (seen[s]) continue;
        n++;
        const stack = [s]; seen[s] = true;
        while (stack.length) { const u = stack.pop(); nbrs(u).forEach((w) => { if (!seen[w]) { seen[w] = true; stack.push(w); } }); }
      }
      return n;
    }
    function getMolecule() {
      if (!V.length) return { ok: false, error: 'Draw a structure first.' };
      const n = pieces();
      if (n > 1) return { ok: false, error: 'Connect everything into one molecule (you have ' + n + ' separate pieces).' };
      const c = chem();
      if (!c) return { ok: false, error: 'The chemistry module did not load. Reload the page and try again.' };
      const atoms = V.map((v) => ({ el: v.el, charge: v.charge || 0, x: v.x, y: v.y }));
      const stereo = stereoFromDrawing(V, E, maxValence);
      Object.keys(stereo.tags).forEach((k) => { atoms[k].chiral = stereo.tags[k].chiral; atoms[k].chiralNbrs = stereo.tags[k].chiralNbrs; });
      const bonds = E.map((e) => ({ a: e.a, b: e.b, order: e.order }));
      try {
        const r = c.fromGraph(atoms, bonds, 'sketch');
        if (r && r.ok) r.stereo = { set: Object.keys(stereo.tags).length, skipped: stereo.skipped.length };
        return r;
      }
      catch (err) { if (root.console) console.warn('NNSketch: fromGraph failed', err); return { ok: false, error: 'I could not read that drawing. Try redrawing it.' }; }
    }
    function emit(result) {
      const r = result || getMolecule();
      listeners.forEach((cb) => { try { cb(r); } catch (err) { if (root.console) console.warn('NNSketch listener', err); } });
    }
    function changed() { redraw(); emit(); }

    /* ----- element picker ----- */
    let pickerArmed = false; // a pointer click counts only if its press began inside the picker
    function openPicker(i) {
      if (!picker) return;
      pickerFor = i; pickerArmed = false;
      picker.innerHTML = '';
      picker.setAttribute('role', 'dialog');
      picker.setAttribute('aria-label', 'Change atom ' + (i + 1));
      ELEMENTS.forEach((el) => {
        const b = document.createElement('button');
        b.type = 'button'; b.textContent = el; b.dataset.el = el;
        b.setAttribute('aria-pressed', String(V[i].el === el));
        b.setAttribute('aria-label', EL_NAME[el]);
        b.style.color = elInfo(el).color2d;
        picker.appendChild(b);
      });
      const del = document.createElement('button');
      del.type = 'button'; del.textContent = 'Del'; del.dataset.el = 'del'; del.className = 'del';
      del.setAttribute('aria-label', 'Delete atom');
      picker.appendChild(del);
      picker.hidden = false;
      // position near the vertex, clamped inside the canvas wrapper
      const pw = picker.offsetWidth, ph = picker.offsetHeight;
      const off = canvas.offsetLeft || 0, offT = canvas.offsetTop || 0;
      const vx = V[i].x, vy = V[i].y, g = 18;
      // beside the atom on whichever side covers the fewest drawn atoms (then has more room);
      // never over the tapped atom itself
      const clampX = (sx) => Math.max(4, Math.min(sx, W - pw - 4)), clampY = (sy) => Math.max(4, Math.min(sy, H - ph - 4));
      const inside = (v, x0, y0, pad) => v.x >= x0 - pad && v.x <= x0 + pw + pad && v.y >= y0 - pad && v.y <= y0 + ph + pad;
      const cand = [
        [vx + g, vy - ph / 2, W - vx], [vx - g - pw, vy - ph / 2, vx],
        [vx - pw / 2, vy + g, H - vy], [vx - pw / 2, vy - g - ph, vy]
      ].map(([sx, sy, room]) => {
        const x0 = clampX(sx), y0 = clampY(sy);
        const cover = V.filter((v, k) => k !== i && inside(v, x0, y0, 6)).length;
        const self = inside(V[i], x0, y0, 12) ? 1 : 0;
        return { x0, y0, score: self * 1e6 + cover * 1e3 - room };
      }).sort((p, q) => p.score - q.score);
      const spot = [cand[0].x0, cand[0].y0];
      let x = spot[0] + off, y = spot[1] + offT;
      picker.style.left = x + 'px'; picker.style.top = y + 'px';
      const first = picker.querySelector('button[aria-pressed="true"]') || picker.querySelector('button');
      if (first) first.focus({ preventScroll: true });
      redraw();
    }
    function closePicker(refocus) {
      if (!picker || picker.hidden) { pickerFor = -1; return; }
      picker.hidden = true; pickerFor = -1;
      if (refocus) canvas.focus({ preventScroll: true });
      redraw();
    }
    function onPickerClick(e) {
      const b = e.target.closest('button[data-el]');
      if (!b || pickerFor < 0) return;
      // The click that follows the opening tap on touch screens lands on whatever button now
      // sits under the finger (often Delete). Ignore pointer clicks whose press started elsewhere.
      if (e.detail > 0 && !pickerArmed) return;
      pickerArmed = false;
      const i = pickerFor, el = b.dataset.el;
      if (el === 'del') { commitStart(); deleteVertex(i); closePicker(true); changed(); return; }
      if (V[i].el === el) { closePicker(true); return; }
      const need = degree(i);
      if (need > maxValence(el)) {
        closePicker(true);
        emit({ ok: false, error: EL_NAME[el] + ' can make only ' + maxValence(el) + ' bond' + (maxValence(el) === 1 ? '' : 's') + ', and this atom has ' + need + '. Remove a bond first.' });
        return;
      }
      commitStart(); V[i].el = el; closePicker(true); changed();
    }
    function onPickerDown() { pickerArmed = true; }
    function onPickerKey(e) { if (e.key === 'Escape') { e.preventDefault(); closePicker(true); } }
    function onDocDown(e) {
      if (!picker || picker.hidden) return;
      if (picker.contains(e.target) || e.target === canvas) return;
      closePicker(false);
    }
    if (picker) {
      picker.addEventListener('click', onPickerClick);
      picker.addEventListener('keydown', onPickerKey);
      picker.addEventListener('pointerdown', onPickerDown);
      document.addEventListener('pointerdown', onDocDown, true);
    }

    /* ----- editing actions ----- */
    function cycleBond(k) {
      const e = E[k];
      const room = (i, delta) => degree(i) + delta <= maxValence(V[i].el);
      let order = e.order, refused = false;
      for (let step = 0; step < 3; step++) {
        order = order % 3 + 1;
        const delta = order - e.order;
        if (room(e.a, delta) && room(e.b, delta)) break;
        refused = true;
      }
      if (refused) { shake = { edge: k, until: performance.now() + 300 }; animate(); }
      if (order === e.order) { redraw(); return; }
      commitStart(); e.order = order; if (order !== 1) { delete e.st; delete e.from; } changed();
    }
    // Wedge/dash tool on an existing bond: line → wedge/dash (narrow end on the more substituted
    // atom, the usual stereocenter; ties go to the end nearer the tap) → flipped → line.
    function stereoBond(k, p) {
      const e = E[k];
      if (e.order !== 1) { shake = { edge: k, until: performance.now() + 300 }; animate(); emit({ ok: false, error: 'Wedges and dashes go on single bonds.' }); return; }
      commitStart();
      if (e.st !== bondTool) {
        const da = nbrs(e.a).length, db = nbrs(e.b).length;
        let from = da > db ? e.a : db > da ? e.b : (dist(p, V[e.a]) <= dist(p, V[e.b]) ? e.a : e.b);
        e.st = bondTool; e.from = from; e.flipped = false;
      } else if (!e.flipped) { e.from = e.from === e.a ? e.b : e.a; e.flipped = true; }
      else { delete e.st; delete e.from; delete e.flipped; }
      changed();
    }
    function tapBond(k, p) { if (bondTool !== 'plain') stereoBond(k, p); else cycleBond(k); }
    function tapFree(p, type) {
      const vi = nearestVertex(p, snapR(type));
      if (vi >= 0) { openPicker(vi); return; }
      const ei = nearestEdge(p, bondR(type));
      if (ei >= 0) { tapBond(ei, p); return; }
      closePicker(false);
    }
    function tapChain(p, type) {
      const vi = nearestVertex(p, snapR(type));
      if (vi >= 0) {
        if (active === vi) { openPicker(vi); return; }
        commitStart();
        if (active != null && active !== vi) addEdge(active, vi);
        active = vi; changed(); return;
      }
      const ei = nearestEdge(p, bondR(type));
      if (ei >= 0) { tapBond(ei, p); return; }
      commitStart();
      let q = p;
      if (active != null) {
        // fixed bond length, direction snapped to multiples of 30 degrees
        const a = V[active];
        const ang = Math.round(Math.atan2(p.y - a.y, p.x - a.x) / (Math.PI / 6)) * (Math.PI / 6);
        q = { x: a.x + bondLength * Math.cos(ang), y: a.y + bondLength * Math.sin(ang) };
        q.x = Math.max(8, Math.min(W - 8, q.x)); q.y = Math.max(8, Math.min(H - 8, q.y));
      }
      let ni = nearestVertex(q, snapR(type));
      if (ni < 0) ni = addVertex(q);
      if (active != null) addEdge(active, ni);
      active = ni;
      changed();
      animate();
    }
    function commitStroke(raw, type) {
      const pts = cleanStroke(raw);
      if (pts.length < 2) { redraw(); return; }
      const R = snapR(type);
      // ring closure: ends near its own start
      let closeRing = false;
      if (pts.length >= 4 && dist(pts[0], pts[pts.length - 1]) < R) { pts.pop(); closeRing = true; }
      commitStart();
      const ids = [];
      pts.forEach((p) => {
        let vi = nearestVertex(p, R);
        if (vi < 0) vi = addVertex(p);
        if (ids[ids.length - 1] !== vi) ids.push(vi);
      });
      if (closeRing && ids.length > 2) ids.push(ids[0]);
      let added = 0;
      const tool = bondTool;
      for (let k = 1; k < ids.length; k++) {
        if (ids[k] !== ids[k - 1] && edgeIndex(ids[k - 1], ids[k]) < 0) {
          bondTool = k === 1 ? tool : 'plain'; // a wedge/dash stroke: only its first bond, narrow end where the stroke began
          addEdge(ids[k - 1], ids[k]); added++;
        }
      }
      bondTool = tool;
      if (!added && ids.length < 2) { history.pop(); redraw(); return; }
      changed();
    }

    /* ----- pointer input ----- */
    let ptr = null;          // the one pointer we are tracking
    let lastPen = -1e9;
    let penDown = false;
    let lpTimer = 0;
    function local(e) { const r = canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; }
    function down(e) {
      if (e.pointerType === 'mouse' && e.button !== 0) return;
      const now = performance.now();
      if (e.pointerType === 'pen') { penDown = true; lastPen = now; }
      // palm rejection: ignore touches while a pen is (or was just) in use
      if (e.pointerType === 'touch' && (penDown || now - lastPen < 1000)) return;
      if (ptr) { // a second finger: cancel the stroke instead of drawing garbage
        if (e.pointerType === 'touch' && ptr.type === 'touch') { ptr = null; stroke = null; clearTimeout(lpTimer); redraw(); }
        return;
      }
      e.preventDefault();
      try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      if (pickerFor >= 0) closePicker(false);
      const p = local(e);
      ptr = { id: e.pointerId, type: e.pointerType, t0: now, start: p, len: 0, last: p, long: false };
      stroke = mode === 'free' ? [p] : null;
      clearTimeout(lpTimer);
      lpTimer = setTimeout(() => {
        if (!ptr || ptr.len > 8) return;
        const vi = nearestVertex(ptr.start, snapR(ptr.type));
        if (vi >= 0) { ptr.long = true; stroke = null; openPicker(vi); }
      }, 500);
      redraw();
    }
    function move(e) {
      if (!ptr || e.pointerId !== ptr.id) {
        if (!ptr && e.pointerType === 'mouse') { // hover feedback
          const h = nearestVertex(local(e), 16);
          const he = h < 0 ? nearestEdge(local(e), 10) : -1;
          canvas.style.cursor = h >= 0 || he >= 0 ? 'pointer' : 'crosshair';
          if (h !== hover) { hover = h; redraw(); }
        }
        return;
      }
      if (e.pointerType === 'pen') lastPen = performance.now();
      const evs = (e.getCoalescedEvents && e.getCoalescedEvents().length) ? e.getCoalescedEvents() : [e];
      evs.forEach((ev) => {
        const p = local(ev);
        ptr.len += dist(p, ptr.last); ptr.last = p;
        if (stroke) stroke.push(p);
      });
      if (ptr.len > 8) clearTimeout(lpTimer);
      if (stroke) redraw();
    }
    function up(e) {
      if (e.pointerType === 'pen') { penDown = false; lastPen = performance.now(); }
      if (!ptr || e.pointerId !== ptr.id) return;
      clearTimeout(lpTimer);
      const p = ptr, raw = stroke;
      ptr = null; stroke = null;
      if (e.type === 'pointercancel' || p.long) { redraw(); return; }
      const tap = p.len < 8 && performance.now() - p.t0 < 400;
      if (mode === 'chain') {
        // In chain mode every short gesture places a point; a slow tap still counts.
        if (p.len < 24) tapChain(p.start, p.type); else redraw();
        return;
      }
      if (tap) { tapFree(p.start, p.type); return; }
      if (raw && raw.length > 1) commitStroke(raw, p.type); else redraw();
    }
    function key(e) {
      if ((e.ctrlKey || e.metaKey) && (e.key === 'z' || e.key === 'Z')) { e.preventDefault(); undo(); }
      else if (e.key === 'Escape') closePicker(true);
    }
    canvas.addEventListener('pointerdown', down);
    canvas.addEventListener('pointermove', move);
    canvas.addEventListener('pointerup', up);
    canvas.addEventListener('pointercancel', up);
    canvas.addEventListener('pointerleave', () => { if (hover >= 0) { hover = -1; redraw(); } });
    canvas.addEventListener('keydown', key);

    /* ----- rendering ----- */
    let raf = 0;
    function animate() { if (!raf) raf = requestAnimationFrame(tick); }
    function tick() {
      raf = 0;
      redraw();
      const now = performance.now();
      if ((shake && now < shake.until) || (mode === 'chain' && active != null && !reduced())) raf = requestAnimationFrame(tick);
      else if (shake && now >= shake.until) { shake = null; redraw(); }
    }
    const reduced = () => !!(root.matchMedia && root.matchMedia('(prefers-reduced-motion: reduce)').matches);

    function hLabel(i) {
      const v = V[i];
      const h = Math.max(0, maxValence(v.el) - degree(i));
      return v.el + (h ? 'H' + (h > 1 ? h : '') : '');
    }
    // Draws e.g. 'NH2' with the element symbol centred on the vertex and the H part
    // trailing to the right; digits are subscripts. A night-coloured knockout sits behind.
    function drawLabel(text, el, x, y, color) {
      const parts = text.match(/\d+|\D+/g) || [];
      let w = 0;
      parts.forEach((t) => { ctx.font = '700 ' + (/\d/.test(t) ? 11 : 15) + 'px ' + FONT; w += ctx.measureText(t).width; });
      ctx.font = '700 15px ' + FONT;
      const elW = ctx.measureText(el).width;
      const left = x - elW / 2;
      ctx.fillStyle = COL.bg;
      roundRect(left - 4, y - 11, w + 8, 22, 9); ctx.fill();
      let cx = left;
      ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; ctx.fillStyle = color;
      parts.forEach((t) => {
        const sub = /\d/.test(t);
        ctx.font = '700 ' + (sub ? 11 : 15) + 'px ' + FONT;
        ctx.fillText(t, cx, y + (sub ? 4 : 0));
        cx += ctx.measureText(t).width;
      });
    }
    function roundRect(x, y, w, h, r) {
      ctx.beginPath();
      ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
      ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
    }
    function valenceBad(i) { return degree(i) > maxValence(V[i].el); }

    function redraw() {
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.fillStyle = COL.bg; ctx.fillRect(0, 0, W, H);
      ctx.fillStyle = COL.grid;
      for (let x = 11; x < W; x += 22) for (let y = 11; y < H; y += 22) ctx.fillRect(x - 1, y - 1, 2, 2);

      if (!V.length && !stroke) {
        // faint dashed ghost zigzag (gone on the first stroke)
        if (mode !== 'chain') {
          const n = 5, sx = 44 * Math.cos(Math.PI / 6), sy = 22, x0 = W / 2 - 2 * sx, y0 = H / 2 + 20;
          ctx.save(); ctx.setLineDash([6, 7]); ctx.strokeStyle = 'rgba(251,243,226,.22)'; ctx.lineWidth = 3; ctx.lineCap = 'round';
          ctx.beginPath();
          for (let k = 0; k < n; k++) { const x = x0 + k * sx, y = y0 + (k % 2 ? -sy : 0); if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
          ctx.stroke(); ctx.restore();
        }
        ctx.fillStyle = COL.dim; ctx.font = '500 14px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.fillText(mode === 'chain' ? 'Click or tap to place atoms' : 'Draw a zigzag here', W / 2, H / 2 - 30);
      }

      const labeled = V.map((v, i) => v.el !== 'C' || (V.length === 1 && !E.length && i === 0));
      const now = performance.now();
      ctx.lineCap = 'round';
      E.forEach((e, k) => {
        let a = V[e.a], b = V[e.b];
        const dx = b.x - a.x, dy = b.y - a.y, L = Math.hypot(dx, dy) || 1, ux = dx / L, uy = dy / L;
        // stop short of heteroatom labels
        const ta = labeled[e.a] ? 11 : 0, tb = labeled[e.b] ? 11 : 0;
        a = { x: a.x + ux * ta, y: a.y + uy * ta }; b = { x: b.x - ux * tb, y: b.y - uy * tb };
        let jitter = 0, color = COL.bond;
        if (shake && shake.edge === k && now < shake.until) { color = COL.err; jitter = Math.sin(now / 25) * 2; }
        const nx = -uy, ny = ux;
        ctx.strokeStyle = color; ctx.lineWidth = 3;
        const line = (p, q, off, shorten) => {
          const sx = (q.x - p.x) * shorten, sy = (q.y - p.y) * shorten;
          ctx.beginPath();
          ctx.moveTo(p.x + sx + nx * (off + jitter), p.y + sy + ny * (off + jitter));
          ctx.lineTo(q.x - sx + nx * (off + jitter), q.y - sy + ny * (off + jitter));
          ctx.stroke();
        };
        if (e.order === 1 && e.st) {
          // narrow end at e.from; wide end 7 px half-width
          const n0 = e.from === e.a ? a : b, w0 = e.from === e.a ? b : a;
          const wx = -(w0.y - n0.y) / L, wy = (w0.x - n0.x) / L, hw = 7;
          ctx.fillStyle = color;
          if (e.st === 'wedge') {
            ctx.beginPath(); ctx.moveTo(n0.x + jitter, n0.y);
            ctx.lineTo(w0.x + wx * hw, w0.y + wy * hw); ctx.lineTo(w0.x - wx * hw, w0.y - wy * hw);
            ctx.closePath(); ctx.fill();
          } else {
            const steps = Math.max(4, Math.round(Math.hypot(w0.x - n0.x, w0.y - n0.y) / 6));
            ctx.lineWidth = 2.2;
            for (let s = 1; s <= steps; s++) {
              const t = s / steps, cx = n0.x + (w0.x - n0.x) * t, cy = n0.y + (w0.y - n0.y) * t, r = 1 + (hw - 1) * t;
              ctx.beginPath(); ctx.moveTo(cx + wx * r + jitter, cy + wy * r); ctx.lineTo(cx - wx * r + jitter, cy - wy * r); ctx.stroke();
            }
            ctx.lineWidth = 3;
          }
        } else if (e.order === 1) line(a, b, 0, 0);
        else if (e.order === 3) { line(a, b, 0, 0); line(a, b, 5.5, 0.1); line(a, b, -5.5, 0.1); }
        else {
          // second line on the side with more neighbours; centred pair when terminal
          const others = nbrs(e.a).filter((j) => j !== e.b).concat(nbrs(e.b).filter((j) => j !== e.a));
          if (!others.length) { line(a, b, 3, 0); line(a, b, -3, 0); }
          else {
            let side = 0;
            others.forEach((j) => { side += (V[j].x - V[e.a].x) * nx + (V[j].y - V[e.a].y) * ny; });
            line(a, b, 0, 0); line(a, b, side >= 0 ? 6 : -6, 0.15);
          }
        }
      });

      V.forEach((v, i) => {
        if (labeled[i]) {
          const text = (V.length === 1 && v.el === 'C') ? 'CH4' : hLabel(i);
          drawLabel(text, v.el, v.x, v.y, elInfo(v.el).color2d);
        } else {
          ctx.fillStyle = COL.mari; ctx.beginPath(); ctx.arc(v.x, v.y, 5, 0, Math.PI * 2); ctx.fill();
          ctx.strokeStyle = COL.bg; ctx.lineWidth = 1.5; ctx.stroke();
        }
        if (valenceBad(i)) { ctx.strokeStyle = COL.err; ctx.lineWidth = 2.5; ctx.beginPath(); ctx.arc(v.x, v.y, 13, 0, Math.PI * 2); ctx.stroke(); }
        if (i === hover || i === pickerFor) { ctx.strokeStyle = 'rgba(255,179,0,.6)'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(v.x, v.y, 11, 0, Math.PI * 2); ctx.stroke(); }
      });

      if (mode === 'chain' && active != null && V[active]) {
        const v = V[active];
        const pulse = reduced() ? 0 : (Math.sin(now / 250) + 1) / 2;
        ctx.strokeStyle = COL.teal; ctx.lineWidth = 2.5;
        ctx.globalAlpha = 0.55 + 0.45 * pulse;
        ctx.beginPath(); ctx.arc(v.x, v.y, 10 + 3 * pulse, 0, Math.PI * 2); ctx.stroke();
        ctx.globalAlpha = 1;
      }

      if (stroke && stroke.length > 1) {
        ctx.strokeStyle = COL.stroke; ctx.lineWidth = 2.5; ctx.lineJoin = 'round';
        ctx.beginPath(); ctx.moveTo(stroke[0].x, stroke[0].y);
        for (let k = 1; k < stroke.length; k++) ctx.lineTo(stroke[k].x, stroke[k].y);
        ctx.stroke();
      }
    }

    /* ----- public handle ----- */
    function undo() {
      if (!history.length) return;
      const s = history.pop();
      V = s.V; E = s.E; active = mode === 'chain' ? (s.active != null && s.active < V.length ? s.active : null) : null;
      closePicker(false); changed();
    }
    function clear() {
      if (!V.length) return;
      commitStart(); V = []; E = []; active = null; closePicker(false); changed();
    }
    function setMode(m) {
      mode = m === 'chain' ? 'chain' : 'free';
      active = null; closePicker(false);
      canvas.style.cursor = 'crosshair';
      redraw();
    }
    function destroy() {
      if (ro) ro.disconnect();
      canvas.removeEventListener('pointerdown', down);
      canvas.removeEventListener('pointermove', move);
      canvas.removeEventListener('pointerup', up);
      canvas.removeEventListener('pointercancel', up);
      canvas.removeEventListener('keydown', key);
      if (picker) {
        picker.removeEventListener('click', onPickerClick);
        picker.removeEventListener('keydown', onPickerKey);
        picker.removeEventListener('pointerdown', onPickerDown);
        document.removeEventListener('pointerdown', onDocDown, true);
      }
      if (raf) cancelAnimationFrame(raf);
    }

    canvas.style.cursor = 'crosshair';
    resize();
    return {
      getMolecule, setMode, undo, clear, redraw, resize, destroy,
      setBondTool(t) { bondTool = t === 'wedge' || t === 'dash' ? t : 'plain'; },
      getBondTool: () => bondTool,
      onChange(cb) { if (typeof cb === 'function') listeners.push(cb); },
      canUndo: () => history.length > 0,
      isEmpty: () => V.length === 0,
      // extra (not in the contract): used by the page's smoke test
      _graph: () => ({ vertices: V.map((v) => Object.assign({}, v)), edges: E.map((e) => Object.assign({}, e)) })
    };
  }

  const api = { attach, _clean: cleanStroke, _stereo: stereoFromDrawing };
  root.NNSketch = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

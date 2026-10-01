/*
 * newman/chain2d.js  —  window.NNChain2D (Newman Navigator, main-chain tracing pad)
 *
 * Draws a typed, preset or drawn molecule as a 2D skeletal structure and reports taps on
 * its atoms. It never decides chemistry: app.js owns the rules and sends a PadState.
 * Contract: newman/CHAIN_SPEC.md section 4. No dependencies, no random numbers.
 *
 * Resolved details:
 *  - The backbone ("as written": heavy atom 0, then the highest-index unvisited carbon
 *    neighbour each step) stops when it enters a ring; the ring is drawn as a regular
 *    polygon and anything beyond it leaves radially and zigzags on as a branch.
 *  - An sp atom (triple bond, or two double bonds) keeps its line straight instead of
 *    zigzagging, the way textbooks draw alkynes.
 *  - A branch whose zigzag direction is not fixed by its parent takes the side with more
 *    room (largest distance to atoms already placed).
 *  - Fused or spiro rings (outside the naming scope) are placed on the shared edge or
 *    atom; anything left unplaced falls back to the branch rule. Good enough to trace.
 */
(function (root) {
  'use strict';

  const COL = {
    bg: '#17131d', grid: '#2d2436', bond: '#fbf3e2', dot: '#9d9386', err: '#e53935', teal: '#14b8a6',
    mari: '#ffb300', ink: '#0f0b13', pink: '#ff5c93', dim: '#9d9386'
  };
  const FONT = 'Inter, -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif';
  const D2R = Math.PI / 180;
  const LABEL_COLOR = { N: '#6f9bff', O: '#ff6b61', F: '#9be15d', Cl: '#2ec27e', Br: '#e0875f', I: '#b784f0', S: '#f2c200', B: '#ffa07a', P: '#ff8c00', H: '#fbf3e2', C: '#ffb300' };

  /* ------------------------------------------------------------------ */
  /* Layout                                                              */
  /* ------------------------------------------------------------------ */
  function heavyGraph(mol) {
    const H = [];
    mol.atoms.forEach((a, i) => { if (a.el !== 'H') H.push(i); });
    const adj = {};
    H.forEach((i) => { adj[i] = []; });
    mol.bonds.forEach((b) => {
      if (adj[b.a] && adj[b.b]) { adj[b.a].push({ n: b.b, o: b.order }); adj[b.b].push({ n: b.a, o: b.order }); }
    });
    H.forEach((i) => adj[i].sort((x, y) => x.n - y.n));
    return { H, adj };
  }
  function hTotal(mol, i) {
    const a = mol.atoms[i];
    if (a.hTotal != null) return a.hTotal;
    let n = a.hCount | 0;
    mol.bonds.forEach((b) => { if ((b.a === i && mol.atoms[b.b].el === 'H') || (b.b === i && mol.atoms[b.a].el === 'H')) n++; });
    return n;
  }

  function layout(mol, opts) {
    opts = opts || {};
    const { H, adj } = heavyGraph(mol);
    const src = opts.src;
    if (src && src.source === 'sketch' && H.every((i) => src.atoms[i] && src.atoms[i].x != null && src.atoms[i].y != null)) {
      const lens = [];
      src.bonds.forEach((b) => {
        const p = src.atoms[b.a], q = src.atoms[b.b];
        if (p.el === 'H' || q.el === 'H') return;
        lens.push(Math.hypot(p.x - q.x, p.y - q.y));
      });
      lens.sort((a, b) => a - b);
      const med = lens.length ? (lens.length % 2 ? lens[(lens.length - 1) / 2] : (lens[lens.length / 2 - 1] + lens[lens.length / 2]) / 2) : 1;
      const s = med > 1e-9 ? 1 / med : 1;
      const pos = {};
      H.forEach((i) => { pos[i] = [src.atoms[i].x * s, src.atoms[i].y * s]; });
      return { pos, L: 1 };
    }
    return { pos: generate(mol, H, adj), L: 1 };
  }

  function generate(mol, H, adj) {
    const pos = {}, parent = {}, dirIn = {}, turn = {}, order = [];
    const isC = (i) => mol.atoms[i].el === 'C';
    const rings = (mol.rings || []).filter((r) => r.every((i) => adj[i]));
    const ringsOf = (i) => rings.filter((r) => r.indexOf(i) >= 0);
    const inRingBond = (a, b) => rings.some((r) => { const i = r.indexOf(a), j = r.indexOf(b); return i >= 0 && j >= 0 && (Math.abs(i - j) === 1 || Math.abs(i - j) === r.length - 1); });
    const linear = (i) => { let d = 0, t = 0; adj[i].forEach((q) => { if (q.o === 2) d++; if (q.o === 3) t++; }); return t > 0 || d >= 2; };
    const ringPlaced = new Set(), ringEntry = new Set();
    const put = (i, p, par, d) => { pos[i] = p; parent[i] = par; dirIn[i] = d; order.push(i); };
    const at = (p, ang, len) => [p[0] + Math.cos(ang * D2R) * (len || 1), p[1] + Math.sin(ang * D2R) * (len || 1)];
    const room = (p, skip) => {
      let m = Infinity;
      Object.keys(pos).forEach((k) => { if (+k === skip) return; const q = pos[k]; m = Math.min(m, Math.hypot(p[0] - q[0], p[1] - q[1])); });
      return m;
    };
    if (!H.length) return pos;

    // place a ring as a polygon. entry atom x already placed, centre in direction `ang` from it
    function placeRing(r, x, ang) {
      const n = r.length, Rc = 1 / (2 * Math.sin(Math.PI / n));
      const c = at(pos[x], ang, Rc);
      const k0 = r.indexOf(x);
      // walk direction: towards the lower-index ring neighbour first (deterministic)
      const nA = r[(k0 + 1) % n], nB = r[(k0 - 1 + n) % n];
      const dir = nA <= nB ? 1 : -1;
      const a0 = ang + 180;
      for (let k = 1; k < n; k++) {
        const ai = r[((k0 + dir * k) % n + n) % n];
        if (pos[ai]) continue;
        const a = a0 + k * 360 / n;
        put(ai, at(c, a, Rc), r[((k0 + dir * (k - 1)) % n + n) % n], a);
        dirIn[ai] = a; // radial, used for exocyclic children
      }
      ringPlaced.add(r);
      dirIn[x] = a0; // children of the entry atom also leave radially
      if (parent[x] != null && parent[x] >= 0) ringEntry.add(x);
      ringCentre.set(r, c);
      fused(r);
    }
    const ringCentre = new Map();
    function fused(r0) {
      rings.forEach((r) => {
        if (ringPlaced.has(r)) return;
        const shared = r.filter((i) => pos[i]);
        if (!shared.length) return;
        const n = r.length, Rc = 1 / (2 * Math.sin(Math.PI / n));
        if (shared.length === 1) {
          const s = shared[0], c0 = ringCentre.get(r0) || [0, 0];
          const ang = Math.atan2(pos[s][1] - c0[1], pos[s][0] - c0[0]) / D2R;
          placeRing(r, s, ang);
          return;
        }
        // shared edge u-v (adjacent in r): centre on the far side from r0's centre
        let u = -1, v = -1;
        for (let k = 0; k < n; k++) { const a = r[k], b = r[(k + 1) % n]; if (pos[a] && pos[b]) { u = a; v = b; break; } }
        if (u < 0) return;
        const c0 = ringCentre.get(r0) || [0, 0];
        const mx = (pos[u][0] + pos[v][0]) / 2, my = (pos[u][1] + pos[v][1]) / 2;
        let nx = -(pos[v][1] - pos[u][1]), ny = pos[v][0] - pos[u][0];
        const L = Math.hypot(nx, ny) || 1; nx /= L; ny /= L;
        if ((mx - c0[0]) * nx + (my - c0[1]) * ny < 0) { nx = -nx; ny = -ny; }
        const ap = 1 / (2 * Math.tan(Math.PI / n));
        const c = [mx + nx * ap, my + ny * ap];
        const ku = r.indexOf(u), kv = r.indexOf(v);
        const dir = (kv - ku + n) % n === 1 ? 1 : -1;
        const av = Math.atan2(pos[v][1] - c[1], pos[v][0] - c[0]) / D2R;
        const au = Math.atan2(pos[u][1] - c[1], pos[u][0] - c[0]) / D2R;
        let step = 360 / n; if (((av - au + 540) % 360) - 180 < 0) step = -step;
        for (let k = 1; k < n - 1; k++) {
          const ai = r[((kv + dir * k) % n + n) % n];
          if (pos[ai]) continue;
          const a = av + k * step;
          put(ai, at(c, a, Rc), r[((kv + dir * (k - 1)) % n + n) % n], a);
        }
        ringPlaced.add(r); ringCentre.set(r, c);
        fused(r);
      });
    }

    // start atom and backbone
    const cs = H.filter(isC);
    const cdeg = (i) => adj[i].filter((q) => isC(q.n)).length;
    let start;
    if (cs.length) {
      const term = cs.filter((i) => cdeg(i) <= 1);
      start = (isC(H[0]) && cdeg(H[0]) <= 1) ? H[0] : (term.length ? term[0] : cs[0]);
    } else start = H[0];
    const backbone = [start];
    if (!ringsOf(start).length) {
      let cur = start;
      const vis = new Set([start]);
      for (;;) {
        const nx = adj[cur].filter((q) => isC(q.n) && !vis.has(q.n)).map((q) => q.n);
        if (!nx.length) break;
        const nxt = nx[nx.length - 1];
        backbone.push(nxt); vis.add(nxt); cur = nxt;
        if (ringsOf(nxt).length) break;
      }
    }
    const contOf = {};
    put(start, [0, 0], -1, null);
    if (ringsOf(start).length) {
      // ring first: polygon around the origin, start atom at the top
      const r = ringsOf(start)[0];
      pos[start] = [0, 0];
      placeRing(r, start, 90);
      dirIn[start] = -90;
    } else {
      // zigzag: -30° (up-right), +30°, -30°, ... ; an sp atom keeps the line straight
      let d = -30;
      for (let k = 1; k < backbone.length; k++) {
        const a = backbone[k - 1], b = backbone[k];
        if (k > 1) d = linear(a) ? dirIn[a] : -dirIn[a];
        contOf[a] = d;
        put(b, at(pos[a], d), a, d);
        if (ringsOf(b).length) { placeRing(ringsOf(b)[0], b, d); }
      }
    }

    // breadth-first: place every remaining neighbour of each placed atom
    const queue = order.slice();
    for (let qi = 0; qi < queue.length; qi++) {
      const i = queue[qi];
      const kids = adj[i].filter((q) => !pos[q.n]).map((q) => q.n);
      if (!kids.length) continue;
      const d = dirIn[i];
      let slots;
      const inRing = ringsOf(i).length > 0;
      if (d == null) {
        // start atom: the backbone leaves at -30°, the rest continue the zigzag backwards
        if (contOf[i] != null) {
          const o = contOf[i];
          slots = kids.length === 1 ? [o - 120] : kids.length === 2 ? [o - 120, o + 120] : [o + 90, o + 180, o - 90];
        } else slots = kids.map((_, k) => -30 + k * 360 / kids.length);
      } else if (inRing) {
        // exocyclic substituents point radially outward; beside the bond we came in on if any
        if (ringEntry.has(i)) slots = [d + 60, d - 60];
        else slots = kids.length === 1 ? [d] : kids.length === 2 ? [d - 35, d + 35] : [d - 50, d, d + 50];
      } else if (contOf[i] != null) {
        const c = contOf[i];
        const f = 2 * d - c;
        slots = kids.length === 1 ? [f] : [f, f + 180, d];
      } else if (linear(i)) {
        slots = [d];
      } else if (kids.length === 1) {
        let s = turn[parent[i]] != null ? -turn[parent[i]] : 0;
        if (!s) {
          const pa = at(pos[i], d + 60), pb = at(pos[i], d - 60);
          s = room(pa, i) >= room(pb, i) - 1e-9 ? 60 : -60;
        }
        turn[i] = s;
        slots = [d + s];
      } else if (kids.length === 2) {
        slots = [d - 60, d + 60];
      } else slots = [d - 90, d, d + 90];
      // kids that are ring atoms enter their ring; others are single atoms
      kids.forEach((x, k) => {
        if (pos[x]) return;
        const ang = slots[Math.min(k, slots.length - 1)] + (k >= slots.length ? 30 * (k - slots.length + 1) : 0);
        put(x, at(pos[i], ang), i, ang);
        queue.push(x);
        const rr = ringsOf(x).filter((r) => !ringPlaced.has(r) && r.indexOf(i) < 0);
        if (rr.length && !inRingBond(i, x)) {
          const before = order.length;
          placeRing(rr[0], x, ang);
          for (let m = before; m < order.length; m++) queue.push(order[m]);
        }
      });
    }
    collisionPass(mol, H, adj, pos, parent, order, inRingBond);
    // tidy rounding so layouts are stable to compare
    H.forEach((i) => { if (pos[i]) pos[i] = [Math.round(pos[i][0] * 1e9) / 1e9, Math.round(pos[i][1] * 1e9) / 1e9]; });
    return pos;
  }

  function collisionPass(mol, H, adj, pos, parent, order, inRingBond) {
    const rank = {}; order.forEach((i, k) => { rank[i] = k; });
    const bonded = (a, b) => adj[a].some((q) => q.n === b);
    const pairs = () => {
      const out = [];
      for (let x = 0; x < H.length; x++) for (let y = x + 1; y < H.length; y++) {
        const a = H[x], b = H[y];
        if (bonded(a, b)) continue;
        if (Math.hypot(pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]) < 0.6) out.push([a, b]);
      }
      return out;
    };
    const side = (x, p) => {
      const seen = new Set([x]), q = [x];
      while (q.length) { const y = q.shift(); for (const { n } of adj[y]) if (n !== p && !seen.has(n)) { seen.add(n); q.push(n); } }
      return seen;
    };
    for (let it = 0; it < 24; it++) {
      const bad = pairs();
      if (!bad.length) return;
      const [a0, b0] = bad[0];
      const late = rank[a0] > rank[b0] ? a0 : b0, other = late === a0 ? b0 : a0;
      // candidate pivots: ancestors of the later atom whose bond to their parent is acyclic
      const cands = [];
      let y = late;
      while (y != null && parent[y] != null && parent[y] >= 0) {
        const p = parent[y];
        if (!inRingBond(p, y)) {
          const S = side(y, p);
          if (!S.has(p) && !S.has(other)) cands.push({ S, p });
        }
        y = p;
      }
      cands.sort((u, v) => u.S.size - v.S.size);
      let best = null;
      for (const c of cands) {
        for (const ang of [60, -60, 120, -120]) {
          const save = {};
          c.S.forEach((i) => { save[i] = pos[i]; });
          const cs = Math.cos(ang * D2R), sn = Math.sin(ang * D2R), o = pos[c.p];
          c.S.forEach((i) => { const dx = pos[i][0] - o[0], dy = pos[i][1] - o[1]; pos[i] = [o[0] + dx * cs - dy * sn, o[1] + dx * sn + dy * cs]; });
          const n = pairs().length;
          if (!best || n < best.n) best = { n, c, ang };
          c.S.forEach((i) => { pos[i] = save[i]; });
          if (n === 0) break;
        }
        if (best && best.n === 0) break;
      }
      if (!best || best.n >= bad.length) return;
      const cs = Math.cos(best.ang * D2R), sn = Math.sin(best.ang * D2R), o = pos[best.c.p];
      best.c.S.forEach((i) => { const dx = pos[i][0] - o[0], dy = pos[i][1] - o[1]; pos[i] = [o[0] + dx * cs - dy * sn, o[1] + dx * sn + dy * cs]; });
    }
  }

  /* ------------------------------------------------------------------ */
  /* Pad: drawing + input                                                 */
  /* ------------------------------------------------------------------ */
  function attach(canvas, h) {
    h = h || {};
    const ctx = canvas.getContext('2d');
    let W = 1, Hh = 1, dpr = 1;
    let mol = null, lay = null, adjG = null, heavy = [];
    let S = emptyState();
    let scr = {}, Lpx = 40;
    let focus = -1, raf = 0, flashT0 = 0, destroyed = false;
    const reduced = () => !!(root.matchMedia && root.matchMedia('(prefers-reduced-motion: reduce)').matches);
    const coarse = () => !!(root.matchMedia && root.matchMedia('(pointer: coarse)').matches);
    const baseH = () => {
      const v = parseFloat((root.getComputedStyle ? root.getComputedStyle(canvas).getPropertyValue('--pad-h') : '') || '');
      return v > 0 ? v : (root.matchMedia && root.matchMedia('(max-width: 860px)').matches ? 220 : 240);
    };

    function emptyState() {
      return { path: [], candidates: [], numbers: {}, numbersDim: {}, branchLabels: {}, ends: [], flash: null, highlight: [], selBond: null, tracing: false };
    }
    function resize() {
      if (destroyed) return;
      dpr = Math.min(root.devicePixelRatio || 1, 2);
      W = canvas.clientWidth || canvas.parentNode && canvas.parentNode.clientWidth || 300;
      // height: the CSS height, grown (never wider) when 30 px bonds would not fit
      let want = baseH();
      if (lay) {
        const bb = bbox();
        const needW = (W - 56) / Math.max(bb.w, 0.5);
        if (needW < 30 || (bb.h * Math.min(56, Math.max(30, needW)) + 56) > want) {
          const L = Math.max(30, Math.min(56, needW));
          want = Math.min(360, Math.max(want, Math.ceil(bb.h * L + 56)));
        }
      }
      canvas.style.height = want + 'px';
      Hh = want;
      canvas.width = Math.round(W * dpr); canvas.height = Math.round(Hh * dpr);
      fit();
      draw();
    }
    function bbox() {
      let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
      heavy.forEach((i) => { const p = lay.pos[i]; if (!p) return; x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); });
      if (!isFinite(x0)) { x0 = y0 = 0; x1 = y1 = 0; }
      return { x0, y0, x1, y1, w: x1 - x0, h: y1 - y0 };
    }
    function fit() {
      if (!lay) return;
      const bb = bbox();
      const fx = (W - 56) / Math.max(bb.w, 0.5), fy = (Hh - 56) / Math.max(bb.h, 0.5);
      Lpx = Math.max(30, Math.min(76, Math.min(fx, fy)));
      const cx = (bb.x0 + bb.x1) / 2, cy = (bb.y0 + bb.y1) / 2;
      scr = {};
      heavy.forEach((i) => { const p = lay.pos[i]; scr[i] = { x: W / 2 + (p[0] - cx) * Lpx, y: Hh / 2 + (p[1] - cy) * Lpx }; });
    }
    const ro = root.ResizeObserver ? new ResizeObserver(() => { const w = canvas.clientWidth; if (Math.abs(w - W) > 0.5) resize(); }) : null;
    if (ro) ro.observe(canvas);

    /* ---- drawing ---- */
    function labelText(i) {
      const a = mol.atoms[i];
      if (a.el === 'C') return null;
      const n = hTotal(mol, i);
      let t = a.el + (n ? 'H' + (n > 1 ? n : '') : '');
      if (a.charge) t += a.charge > 0 ? '+' : '−';
      return t;
    }
    function roundRect(x, y, w, hh, r) {
      ctx.beginPath();
      ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + hh, r); ctx.arcTo(x + w, y + hh, x, y + hh, r);
      ctx.arcTo(x, y + hh, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
    }
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
    function disc(x, y, r, fill, stroke, lw) {
      ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2);
      if (fill) { ctx.fillStyle = fill; ctx.fill(); }
      if (stroke) { ctx.strokeStyle = stroke; ctx.lineWidth = lw || 2; ctx.stroke(); }
    }
    function ringCentreOf(a, b) {
      const rs = (mol.rings || []).filter((r) => r.indexOf(a) >= 0 && r.indexOf(b) >= 0);
      if (!rs.length) return null;
      let x = 0, y = 0; rs[0].forEach((i) => { x += scr[i].x; y += scr[i].y; });
      return { x: x / rs[0].length, y: y / rs[0].length };
    }
    function draw() {
      if (destroyed) return;
      const now = (root.performance && performance.now()) || 0;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.fillStyle = COL.bg; ctx.fillRect(0, 0, W, Hh);
      // a faint grid, so the molecule and its tappable carbons stand out
      ctx.fillStyle = 'rgba(45,36,54,.55)';
      for (let x = 11; x < W; x += 22) for (let y = 11; y < Hh; y += 22) ctx.fillRect(x - 1, y - 1, 2, 2);
      if (!mol || !lay) return;
      const fl = S.flash, flAge = fl ? now - flashT0 : 1e9;
      const flashOn = fl && (fl.kind === 'err' ? flAge < 300 : flAge < 1200);
      let shake = 0;
      if (flashOn && fl.kind === 'err' && !reduced()) shake = Math.sin(flAge / 22) * 4 * (1 - flAge / 300);
      const P = (i) => { const p = scr[i]; return (flashOn && fl.kind === 'err' && fl.atoms.indexOf(i) >= 0) ? { x: p.x + shake, y: p.y } : p; };
      const labeled = {};
      heavy.forEach((i) => { labeled[i] = mol.atoms[i].el !== 'C' || heavy.length === 1; });
      const heavyBonds = mol.bonds.filter((b) => scr[b.a] && scr[b.b]);
      ctx.lineCap = 'round'; ctx.lineJoin = 'round';
      // selected bond (marigold, under everything)
      if (S.selBond && scr[S.selBond[0]] && scr[S.selBond[1]]) {
        const a = P(S.selBond[0]), b = P(S.selBond[1]);
        ctx.strokeStyle = 'rgba(255,179,0,.35)'; ctx.lineWidth = 14;
        ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
      }
      // pink highlight halo
      const hi = new Set(S.highlight || []);
      if (hi.size) {
        ctx.strokeStyle = 'rgba(255,92,147,.55)'; ctx.lineWidth = 12;
        heavyBonds.forEach((b) => { if (hi.has(b.a) && hi.has(b.b)) { const a = P(b.a), c = P(b.b); ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(c.x, c.y); ctx.stroke(); } });
        hi.forEach((i) => { if (scr[i]) { const p = P(i); disc(p.x, p.y, labeled[i] ? 14 : 9, 'rgba(255,92,147,.55)'); } });
      }
      // traced path (teal, under the bonds)
      if (S.path.length > 1) {
        ctx.strokeStyle = COL.teal; ctx.lineWidth = 9;
        ctx.beginPath();
        S.path.forEach((i, k) => { const p = P(i); if (k) ctx.lineTo(p.x, p.y); else ctx.moveTo(p.x, p.y); });
        ctx.stroke();
      }
      // bonds
      heavyBonds.forEach((e) => {
        let a = P(e.a), b = P(e.b);
        const dx = b.x - a.x, dy = b.y - a.y, L = Math.hypot(dx, dy) || 1, ux = dx / L, uy = dy / L;
        const ta = labeled[e.a] ? 12 : 0, tb = labeled[e.b] ? 12 : 0;
        a = { x: a.x + ux * ta, y: a.y + uy * ta }; b = { x: b.x - ux * tb, y: b.y - uy * tb };
        const nx = -uy, ny = ux;
        const isSel = S.selBond && ((S.selBond[0] === e.a && S.selBond[1] === e.b) || (S.selBond[0] === e.b && S.selBond[1] === e.a));
        ctx.strokeStyle = isSel ? COL.mari : COL.bond; ctx.lineWidth = isSel ? 4 : 3;
        const line = (off, shorten) => {
          const sx = (b.x - a.x) * shorten, sy = (b.y - a.y) * shorten;
          ctx.beginPath(); ctx.moveTo(a.x + sx + nx * off, a.y + sy + ny * off); ctx.lineTo(b.x - sx + nx * off, b.y - sy + ny * off); ctx.stroke();
        };
        if (e.order === 1) line(0, 0);
        else if (e.order === 3) { line(0, 0); line(5.5, 0.1); line(-5.5, 0.1); }
        else {
          const rc = ringCentreOf(e.a, e.b);
          let side = 0;
          if (rc) side = (rc.x - a.x) * nx + (rc.y - a.y) * ny;
          else {
            (adjG[e.a] || []).concat(adjG[e.b] || []).forEach((q) => { if (q.n !== e.a && q.n !== e.b && scr[q.n]) side += (scr[q.n].x - scr[e.a].x) * nx + (scr[q.n].y - scr[e.a].y) * ny; });
          }
          if (!rc && Math.abs(side) < 1e-6) { line(3, 0); line(-3, 0); }
          else { line(0, 0); line(side >= 0 ? 6 : -6, 0.15); }
        }
      });
      // atoms
      const inPath = new Set(S.path);
      heavy.forEach((i) => {
        const p = P(i);
        if (labeled[i]) {
          const t = heavy.length === 1 && mol.atoms[i].el === 'C' ? 'CH' + Math.max(0, hTotal(mol, i)) : labelText(i);
          drawLabel(t, mol.atoms[i].el, p.x, p.y, LABEL_COLOR[mol.atoms[i].el] || COL.bond);
        } else if (!inPath.has(i)) disc(p.x, p.y, 6, '#2d2436', COL.bond, 2);   // a visible, tappable carbon
      });
      // traced carbons
      S.path.forEach((i, k) => {
        const p = P(i);
        disc(p.x, p.y, 8, COL.teal, COL.ink, 2);
        if (k === S.path.length - 1 && S.tracing) {
          const t = reduced() ? 0.5 : (now % 1200) / 1200;
          ctx.globalAlpha = reduced() ? 0.8 : 1 - t;
          disc(p.x, p.y, 11 + (reduced() ? 2 : t * 7), null, COL.teal, 2.5);
          ctx.globalAlpha = 1;
        }
      });
      // candidates
      if (S.candidates.length) {
        ctx.save(); ctx.setLineDash([4, 4]);
        S.candidates.forEach((i) => { if (scr[i]) { const p = P(i); disc(p.x, p.y, 14, null, COL.mari, 2); } });
        ctx.restore();
      }
      // numbers
      const num = (map, alpha) => {
        Object.keys(map || {}).forEach((k) => {
          const i = +k; if (!scr[i]) return;
          const p = P(i);
          ctx.globalAlpha = alpha;
          disc(p.x, p.y, 11, COL.mari, COL.ink, 2);
          ctx.fillStyle = COL.ink; ctx.font = '700 12px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          ctx.fillText(String(map[k]), p.x, p.y + 0.5);
          ctx.globalAlpha = 1;
        });
      };
      num(S.numbersDim, 0.45);
      num(S.numbers, 1);
      // branch labels (C5 ...)
      Object.keys(S.branchLabels || {}).forEach((k) => {
        const i = +k; if (!scr[i]) return;
        const p = P(i), off = awayDir(i);
        ctx.fillStyle = COL.dim; ctx.font = '600 11px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.fillText(S.branchLabels[k], p.x + off.x * 16, p.y + off.y * 16);
      });
      // end buttons
      (S.ends || []).forEach((i) => {
        if (!scr[i]) return;
        const p = P(i);
        disc(p.x + 2, p.y + 2, 15, COL.err);
        disc(p.x, p.y, 15, COL.mari, COL.ink, 2);
        ctx.fillStyle = COL.ink; ctx.font = '700 12px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.fillText('1?', p.x, p.y + 0.5);
      });
      // flash
      if (flashOn) {
        if (fl.kind === 'err') fl.atoms.forEach((i) => { if (scr[i]) { const p = P(i); disc(p.x, p.y, 14, null, COL.err, 3); } });
        else {
          const t = (flAge % 400) / 400;
          ctx.globalAlpha = reduced() ? 1 : 1 - t * 0.7;
          fl.atoms.forEach((i) => { if (scr[i]) { const p = P(i); disc(p.x, p.y, 13 + (reduced() ? 2 : t * 6), null, COL.mari, 3); } });
          ctx.globalAlpha = 1;
        }
      }
      // keyboard focus
      if (focus >= 0 && scr[focus] && root.document && document.activeElement === canvas) {
        const p = P(focus);
        ctx.save(); ctx.setLineDash([4, 3]); disc(p.x, p.y, 15, null, COL.mari, 2); ctx.restore();
      }
      const animate = (flashOn) || (S.tracing && S.path.length && !reduced());
      if (animate && !raf) raf = root.requestAnimationFrame(() => { raf = 0; draw(); });
    }
    function awayDir(i) {
      let x = 0, y = 0;
      (adjG[i] || []).forEach((q) => { if (scr[q.n]) { x += scr[i].x - scr[q.n].x; y += scr[i].y - scr[q.n].y; } });
      const L = Math.hypot(x, y);
      return L > 1e-6 ? { x: x / L, y: y / L } : { x: 0, y: -1 };
    }

    /* ---- input ---- */
    function hitR() { return Math.max(coarse() ? 22 : 18, Math.min(30, 0.6 * Lpx)); }
    function local(e) { const r = canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; }
    function nearest(p, rad, list) {
      let best = -1, bd = Infinity;
      (list || heavy).forEach((i) => { const s = scr[i]; if (!s) return; const d = Math.hypot(s.x - p.x, s.y - p.y); if (d < bd) { bd = d; best = i; } });
      return bd <= rad ? best : -1;
    }
    function hit(p) {
      if (S.ends && S.ends.length) {
        const e = nearest(p, Math.max(22, hitR()), S.ends);
        if (e >= 0) return { atom: e, end: true };
      }
      const i = nearest(p, hitR());
      return { atom: i >= 0 ? i : null, end: false };
    }
    let down = null;
    canvas.addEventListener('pointerdown', (e) => {
      if (!mol) return;
      const p = local(e);
      down = { id: e.pointerId, x: p.x, y: p.y, t: (root.performance && performance.now()) || 0, moved: 0, dragged: false, type: e.pointerType };
      if (S.tracing) { try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ } }
    });
    canvas.addEventListener('pointermove', (e) => {
      if (!down || e.pointerId !== down.id) return;
      const p = local(e);
      down.moved = Math.max(down.moved, Math.hypot(p.x - down.x, p.y - down.y));
      if (!S.tracing || down.moved < 10) return;
      // drag tracing: entering a valid next carbon appends it
      const i = nearest(p, hitR() * 0.8);
      if (i < 0) return;
      const last = S.path[S.path.length - 1];
      if (!down.dragged && S.path.length === 0) { down.dragged = true; const s0 = nearest({ x: down.x, y: down.y }, hitR()); if (s0 >= 0 && h.onTap) h.onTap(s0, { x: down.x, y: down.y, drag: true }); return; }
      if (!down.dragged && last != null) {
        const s0 = nearest({ x: down.x, y: down.y }, hitR());
        down.dragged = true;
        if (s0 >= 0 && s0 !== last && h.onTap) h.onTap(s0, { x: down.x, y: down.y, drag: true });
      }
      const cur = S.path[S.path.length - 1];
      if (i !== cur && S.path.indexOf(i) < 0 && (S.candidates || []).indexOf(i) >= 0 && h.onTap) h.onTap(i, { x: p.x, y: p.y, drag: true });
    });
    const up = (e) => {
      if (!down || e.pointerId !== down.id) return;
      const d = down; down = null;
      try { canvas.releasePointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      if (e.type !== 'pointerup' || d.dragged) return;
      const dt = ((root.performance && performance.now()) || 0) - d.t;
      if (d.moved >= 10 || dt >= 500) return;
      const p = local(e), r = hit(p);
      if (r.end && h.onEndTap) { h.onEndTap(r.atom); return; }
      if (h.onTap) h.onTap(r.atom, { x: p.x, y: p.y });
    };
    canvas.addEventListener('pointerup', up);
    canvas.addEventListener('pointercancel', up);
    canvas.addEventListener('keydown', (e) => {
      if (!mol || !heavy.length) return;
      const k = e.key;
      if (k === 'ArrowLeft' || k === 'ArrowRight' || k === 'ArrowUp' || k === 'ArrowDown') {
        e.preventDefault();
        if (focus < 0) { focus = heavy.slice().sort((a, b) => scr[a].x - scr[b].x || scr[a].y - scr[b].y)[0]; draw(); return; }
        const f = scr[focus];
        let best = -1, bd = Infinity;
        heavy.forEach((i) => {
          if (i === focus) return;
          const s = scr[i], dx = s.x - f.x, dy = s.y - f.y;
          let ok, d;
          // nearest atom in that direction, a little extra cost for being off-axis
          const r = Math.hypot(dx, dy);
          if (k === 'ArrowRight') { ok = dx > 1; d = r + Math.abs(dy) * 0.5; }
          else if (k === 'ArrowLeft') { ok = dx < -1; d = r + Math.abs(dy) * 0.5; }
          else if (k === 'ArrowUp') { ok = dy < -1; d = r + Math.abs(dx) * 0.5; }
          else { ok = dy > 1; d = r + Math.abs(dx) * 0.5; }
          if (ok && d < bd) { bd = d; best = i; }
        });
        if (best >= 0) { focus = best; draw(); }
      } else if (k === 'Enter' || k === ' ') {
        if (focus < 0) return;
        e.preventDefault();
        if ((S.ends || []).indexOf(focus) >= 0 && h.onEndTap) h.onEndTap(focus);
        else if (h.onTap) h.onTap(focus, { x: scr[focus].x, y: scr[focus].y, key: true });
      } else if (k === 'Backspace') {
        if (h.onUndo) { e.preventDefault(); h.onUndo(); }
      }
    });
    canvas.addEventListener('focus', () => draw());
    canvas.addEventListener('blur', () => draw());

    return {
      setMol(m, l) {
        mol = m; lay = l || layout(m);
        const g = heavyGraph(m); adjG = g.adj; heavy = g.H.filter((i) => lay.pos[i]);
        focus = -1; S = emptyState();
        resize();
      },
      setState(s) {
        const prevFlash = S.flash;
        S = Object.assign(emptyState(), s || {});
        if (S.flash && S.flash !== prevFlash) flashT0 = (root.performance && performance.now()) || 0;
        canvas.classList.toggle('tracing', !!S.tracing);
        draw();
      },
      resize,
      destroy() { destroyed = true; if (ro) ro.disconnect(); if (raf) root.cancelAnimationFrame(raf); },
      focusAtom(i) { focus = i; draw(); },
      screenPos(i) { return scr[i] ? { x: scr[i].x, y: scr[i].y } : null; },
      get bondPx() { return Lpx; }
    };
  }

  const api = { layout, attach };
  root.NNChain2D = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

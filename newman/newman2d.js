/*
 * newman/newman2d.js  —  window.NNNewman2D (Newman Navigator, 2D drawing)
 *
 * Draws a Newman projection (side panel and the translucent overlay on top of the
 * 3D view) and the torsional-energy plot, in the site's rickshaw style.
 * Contract: newman/SPEC.md section 6. Pure canvas 2D, no dependencies.
 *
 * Angle convention: degrees clockwise from 12 o'clock, as the viewer sees the page.
 * Screen point for angle a at radius R: (cx + R sin a, cy − R cos a).
 *
 * Resolved details (signatures unchanged):
 *  - Overlay style keeps the substituent labels (smaller, with a dark knockout so they
 *    read over the 3D model) but draws no atom labels, as the spec asks.
 *  - A Sub may carry `priority: true`; that is honoured in addition to `highlight`.
 *  - The lock glyph for a fixed bond sits in the bottom-right corner (2.05 r out) so it
 *    never covers a substituent label.
 *  - Front lone pairs are drawn as a dot pair at 0.62 r on their ray (the back ones sit
 *    outside the circle at 1.3 r with a short dim stub, per spec).
 */
(function (root) {
  'use strict';

  const C = {
    night: '#17131d', night2: '#221b2a', line: '#3d3346', cream: '#fbf3e2', muted: '#c9bfae',
    dim: '#9d9386', mari: '#ffb300', pink: '#ff5c93', teal: '#14b8a6', ink: '#0f0b13', red: '#e53935'
  };
  const FONT = 'Inter, -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif';
  const RAD = Math.PI / 180;

  const px = (cx, a, R) => cx + R * Math.sin(a * RAD);
  const py = (cy, a, R) => cy - R * Math.cos(a * RAD);

  /* ---------- chemical text with subscript digits ---------- */
  // Splits 'CH(CH3)2' into runs; digits become subscripts (0.72 size, +3 px baseline).
  function runs(text) {
    const out = [];
    for (const ch of String(text)) {
      const sub = ch >= '0' && ch <= '9';
      const last = out[out.length - 1];
      if (last && last.sub === sub) last.t += ch; else out.push({ t: ch, sub });
    }
    return out;
  }
  function chemWidth(ctx, text, size) {
    let w = 0;
    for (const r of runs(text)) {
      ctx.font = '700 ' + (r.sub ? size * 0.72 : size) + 'px ' + FONT;
      w += ctx.measureText(r.t).width;
    }
    return w;
  }
  // Draws text centred on (x, y) (vertical middle of the main glyphs).
  function chemText(ctx, text, x, y, size, color) {
    const w = chemWidth(ctx, text, size);
    let cx = x - w / 2;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    ctx.fillStyle = color;
    for (const r of runs(text)) {
      const s = r.sub ? size * 0.72 : size;
      ctx.font = '700 ' + s + 'px ' + FONT;
      ctx.fillText(r.t, cx, y + (r.sub ? size * 0.23 + 1 : 0));
      cx += ctx.measureText(r.t).width;
    }
    return w;
  }
  function roundRect(ctx, x, y, w, h, r) {
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.arcTo(x + w, y, x + w, y + h, r);
    ctx.arcTo(x + w, y + h, x, y + h, r);
    ctx.arcTo(x, y + h, x, y, r);
    ctx.arcTo(x, y, x + w, y, r);
    ctx.closePath();
  }
  // Label with an optional dark knockout pill and optional pink priority halo. A small
  // grey tag (ax / eq) can follow the text.
  function labelSize(ctx, text, size, tag) {
    let w = chemWidth(ctx, text, size);
    if (tag) { ctx.font = '600 ' + Math.round(size * 0.72) + 'px ' + FONT; w += 3 + ctx.measureText(tag).width; }
    return { w: w + 10, h: size + 6 };
  }
  function label(ctx, text, x, y, size, color, opts) {
    if (!text) return;
    const sz = labelSize(ctx, text, size, opts.tag), w = sz.w - 10, h = sz.h, padX = 5;
    if (opts.knockout || opts.priority) {
      roundRect(ctx, x - w / 2 - padX, y - h / 2, w + 2 * padX, h, h / 2);
      ctx.fillStyle = opts.knockout || 'rgba(23,19,29,.85)';
      ctx.fill();
      if (opts.priority) { ctx.lineWidth = 2.5; ctx.strokeStyle = C.pink; ctx.stroke(); }
    }
    const tw = chemWidth(ctx, text, size);
    chemText(ctx, text, x - w / 2 + tw / 2, y, size, color);
    if (opts.tag) {
      ctx.font = '600 ' + Math.round(size * 0.72) + 'px ' + FONT;
      ctx.fillStyle = C.dim; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText(opts.tag, x - w / 2 + tw + 3, y + 1);
    }
  }
  // Distance from a point to the edge of a w x h box centred on it, along angle a.
  function extent(a, w, h) {
    const sn = Math.abs(Math.sin(a * RAD)), cs = Math.abs(Math.cos(a * RAD));
    return Math.min(sn > 1e-6 ? (w / 2) / sn : 1e9, cs > 1e-6 ? (h / 2) / cs : 1e9);
  }

  /* ---------- lock glyph (drawn with paths, no emoji) ---------- */
  function lock(ctx, x, y, s) {
    ctx.save();
    ctx.strokeStyle = C.muted; ctx.fillStyle = C.muted; ctx.lineWidth = Math.max(1.5, s * 0.14);
    ctx.beginPath(); ctx.arc(x, y - s * 0.15, s * 0.3, Math.PI, 0); ctx.stroke();
    roundRect(ctx, x - s * 0.45, y - s * 0.15, s * 0.9, s * 0.7, s * 0.12); ctx.fill();
    ctx.fillStyle = C.night; ctx.beginPath(); ctx.arc(x, y + s * 0.18, s * 0.1, 0, 7); ctx.fill();
    ctx.restore();
  }

  const wrap180 = (d) => { d = ((d % 360) + 360) % 360; return d > 180 ? d - 360 : d; };

  /* ---------- Newman projection ---------- */
  // Textbook convention: front bonds meet at the centre and run out past the circle; back
  // bonds start at the circle's edge. Both end at 1.5 r, with the label on the same line.
  // When a back bond nearly eclipses a front bond it is drawn 9 degrees off so both show.
  function draw(ctx, d) {
    const { cx, cy } = d;
    const alpha = d.alpha == null ? 1 : d.alpha;
    const r = Math.max(8, d.r * (d.scale == null ? 1 : d.scale));
    const overlay = d.style === 'overlay';
    const lw = 3;
    const size = overlay ? Math.max(11, Math.min(13, r * 0.2)) : 13;
    const L = 1.5 * r;
    const front = d.frontSubs || [], back = d.backSubs || [];
    const hl = d.highlight || [-1, -1];
    const isPri = (list, k, side) => !!(list[k] && (list[k].priority || hl[side] === k));

    // drawn angle for back substituents (offset when eclipsing a front one)
    const ecl = back.map((s) => {
      let best = null, fk = -1;
      front.forEach((f, k) => { const dd = wrap180(s.angle - f.angle); if (!f.lp && Math.abs(dd) < 10 && (!best || Math.abs(dd) < Math.abs(best))) { best = dd; fk = k; } });
      return best == null ? null : { d: best, fk, sgn: best >= 0 ? 1 : -1 };
    });
    const backDraw = back.map((s, k) => (ecl[k] ? front[ecl[k].fk].angle + ecl[k].sgn * 9 : s.angle));

    ctx.save();
    ctx.globalAlpha = alpha;
    ctx.lineCap = 'round';

    const stroke = (x0, y0, x1, y1, color) => {
      ctx.strokeStyle = C.ink; ctx.lineWidth = lw + 3;
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
      ctx.strokeStyle = color || C.cream; ctx.lineWidth = lw;
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    };

    // Back bonds first: from the circle edge out to L.
    back.forEach((s, k) => {
      const a = backDraw[k];
      if (s.lp) {
        ctx.strokeStyle = 'rgba(157,147,134,.55)'; ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.moveTo(px(cx, a, r), py(cy, a, r));
        ctx.lineTo(px(cx, a, r * 1.15), py(cy, a, r * 1.15)); ctx.stroke();
        lpDots(ctx, cx, cy, a, r * 1.3, s.color || C.dim);
        return;
      }
      stroke(px(cx, a, r), py(cy, a, r), px(cx, a, L), py(cy, a, L), C.cream);
    });

    // The back carbon: the circle (a translucent disc over the 3D model).
    ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2);
    ctx.fillStyle = overlay ? 'rgba(23,19,29,.62)' : C.night2;
    ctx.fill();
    ctx.lineWidth = 3; ctx.strokeStyle = C.ink; ctx.stroke();
    ctx.lineWidth = 2.5;
    ctx.strokeStyle = C.cream;
    if (d.rotatable === false) ctx.setLineDash([6, 5]);
    ctx.stroke();
    ctx.setLineDash([]);

    // Dihedral arc between the two priority groups (under the front bonds).
    const fi = hl[0], bi = hl[1];
    let arcLabel = null;
    if (d.dihedral != null && front[fi] && back[bi]) {
      const a0 = front[fi].angle, dih = d.dihedral;
      const ra = r * 0.42;
      if (Math.abs(dih) > 3) {
        ctx.strokeStyle = C.mari; ctx.lineWidth = 2.2;
        ctx.beginPath();
        ctx.arc(cx, cy, ra, (a0 - 90) * RAD, (a0 + dih - 90) * RAD, dih < 0);
        ctx.stroke();
      }
      // angle label on the arc's bisector, inside the circle between front bonds; when the
      // bisector runs along a front bond (eclipsed), use the opposite side
      let mid = a0 + dih / 2;
      if (front.some((f) => !f.lp && Math.abs(wrap180(f.angle - mid)) < 22)) mid += 180;
      let dr = Math.round(dih); if (dr <= -180) dr = 180;
      arcLabel = { a: mid, txt: (dr < 0 ? '−' : '') + Math.abs(dr) + '°' };
    }

    // Front bonds: centre out past the circle to L, over everything else.
    front.forEach((s) => {
      if (s.lp) { lpDots(ctx, cx, cy, s.angle, r * 0.62, s.color || C.dim); return; }
      stroke(cx, cy, px(cx, s.angle, L), py(cy, s.angle, L), C.cream);
    });
    ctx.fillStyle = C.cream;
    ctx.beginPath(); ctx.arc(cx, cy, 4, 0, Math.PI * 2); ctx.fill();
    ctx.lineWidth = 1.5; ctx.strokeStyle = C.ink; ctx.stroke();

    if (arcLabel) {
      ctx.font = '700 12px ' + FONT;
      const tw = ctx.measureText(arcLabel.txt).width;
      const tr = r * 0.7;
      const tx = px(cx, arcLabel.a, tr), ty = py(cy, arcLabel.a, tr);
      roundRect(ctx, tx - tw / 2 - 4, ty - 8, tw + 8, 16, 8);
      ctx.fillStyle = 'rgba(15,11,19,.9)'; ctx.fill();
      ctx.fillStyle = C.mari; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(arcLabel.txt, tx, ty + 0.5);
    }

    // Labels last, on the line of their bond just past its end.
    const ko = overlay ? 'rgba(15,11,19,.86)' : 'rgba(34,27,42,.94)';
    const place = (s, a, pri) => {
      const sz = labelSize(ctx, s.label, size, s.tag);
      const R = L + 4 + extent(a, sz.w, sz.h);
      label(ctx, s.label, px(cx, a, R), py(cy, a, R), size, s.color || C.cream, { knockout: ko, priority: pri, tag: s.tag });
    };
    // eclipsed pairs: spread the two labels apart around the shared direction
    const frontLab = front.map((f) => f.angle), backLab = backDraw.slice();
    ecl.forEach((e, k) => {
      if (!e || !back[k].label || !front[e.fk].label) return;
      const f = front[e.fk];
      const wb = labelSize(ctx, back[k].label, size, back[k].tag).w, wf = labelSize(ctx, f.label, size, f.tag).w;
      const sep = ((wb + wf) / 2 + 4) / (L + 16) / RAD;
      const s1 = Math.max(9, sep / 2 + 2), s2 = Math.max(0, sep - s1);
      backLab[k] = f.angle + e.sgn * s1; frontLab[e.fk] = f.angle - e.sgn * s2;
    });
    back.forEach((s, k) => { if (!s.lp && s.label) place(s, backLab[k], isPri(back, k, 1)); });
    front.forEach((s, k) => { if (!s.lp && s.label) place(s, frontLab[k], isPri(front, k, 0)); });

    if (d.rotatable === false && d.style !== 'overlay') lock(ctx, cx + r * 1.7, cy + r * 1.7, Math.max(12, r * 0.22)); // corner, clear of labels
    ctx.restore();
  }

  function lpDots(ctx, cx, cy, a, R, color) {
    const x = px(cx, a, R), y = py(cy, a, R);
    const tx = Math.cos(a * RAD), ty = Math.sin(a * RAD); // perpendicular to the ray
    ctx.fillStyle = color;
    [-1, 1].forEach((k) => { ctx.beginPath(); ctx.arc(x + k * 4 * tx, y + k * 4 * ty, 2.5, 0, Math.PI * 2); ctx.fill(); });
  }

  /* ---------- torsional energy plot ---------- */
  function drawEnergy(ctx, d) {
    const W = d.w, H = d.h;
    const L = 36, B = 36, T = 30, Rm = 12; // bottom margin holds tick labels plus the axis title
    const pw = Math.max(10, W - L - Rm), ph = Math.max(10, H - T - B);
    const pts = d.points || [];
    let max = 0;
    pts.forEach((p) => { if (p.kJ > max) max = p.kJ; });
    const ymax = Math.ceil(max / 5) * 5 + 5;
    const X = (deg) => L + (deg / 360) * pw;
    const Y = (kJ) => T + ph - (kJ / ymax) * ph;

    ctx.save();
    ctx.clearRect(0, 0, W, H);
    ctx.font = '500 11px ' + FONT;
    ctx.textBaseline = 'middle';

    // grid + axes
    ctx.strokeStyle = C.line; ctx.lineWidth = 1;
    const ystep = ymax > 30 ? 10 : 5;
    for (let v = 0; v <= ymax; v += ystep) {
      ctx.beginPath(); ctx.moveTo(L, Y(v) + 0.5); ctx.lineTo(L + pw, Y(v) + 0.5); ctx.stroke();
      ctx.fillStyle = C.dim; ctx.textAlign = 'right'; ctx.fillText(String(v), L - 5, Y(v));
    }
    for (let deg = 0; deg <= 360; deg += 60) {
      ctx.beginPath(); ctx.moveTo(X(deg) + 0.5, T); ctx.lineTo(X(deg) + 0.5, T + ph); ctx.stroke();
      ctx.fillStyle = C.dim; ctx.textAlign = 'center'; ctx.fillText(deg + '°', X(deg), T + ph + 10);
    }
    // axis titles
    ctx.fillStyle = C.dim; ctx.font = '500 11px ' + FONT;
    ctx.textAlign = 'left'; ctx.fillText(d.yLabel || 'energy (kJ/mol)', 4, 7);
    ctx.textAlign = 'right'; ctx.fillText(d.xLabel || 'dihedral', L + pw, T + ph + 27);

    // curve
    if (pts.length) {
      ctx.strokeStyle = C.teal; ctx.lineWidth = 2.5; ctx.lineJoin = 'round';
      ctx.beginPath();
      pts.forEach((p, k) => { if (k) ctx.lineTo(X(p.deg), Y(p.kJ)); else ctx.moveTo(X(p.deg), Y(p.kJ)); });
      ctx.stroke();
    }

    // Stationary point names, all above their point: maxima just above the peak, minima
    // about 12 px above the curve's minimum. Labels at the plot edges move inward. A label
    // that would collide with another (tight plots on phones) or with the value tag is skipped.
    const cur = d.current ? { x: X(((d.current.deg % 360) + 360) % 360), y: Y(d.current.kJ) } : null;
    ctx.font = '700 11px ' + FONT;
    const tagTxt = d.current ? d.current.kJ.toFixed(1) + ' kJ/mol' : '';
    const tagW = d.current ? ctx.measureText(tagTxt).width : 0;
    ctx.font = '500 11px ' + FONT;
    const boxes = [];
    const hit = (b) => boxes.some((o) => b.x0 < o.x1 + 3 && b.x1 > o.x0 - 3 && b.y0 < o.y1 && b.y1 > o.y0);
    const labels = [];
    (d.stationary || []).slice().sort((a, b) => (a.kind === b.kind ? 0 : a.kind === 'max' ? -1 : 1)).forEach((s) => {
      if (!s.name) return;
      const x = X(s.deg), w = ctx.measureText(s.name).width;
      const xl = Math.min(Math.max(x, L + w / 2 + 2), L + pw - w / 2);
      const y = s.kind === 'max' ? Math.max(T - 10, Y(s.kJ) - 10) : Y(s.kJ) - 12;
      const box = { x0: xl - w / 2, x1: xl + w / 2, y0: y - 7, y1: y + 7 };
      if (hit(box)) return;
      boxes.push(box); labels.push({ t: s.name, x: xl, y });
    });
    // value tag: on the side of the marker with room and no label
    let tag = null;
    if (cur) {
      const ty = Math.max(T - 4, cur.y - 16);
      const opts = [cur.x + 12, cur.x - 12 - tagW - 8].map((x0) => ({ x0, x1: x0 + tagW + 8, y0: ty - 8, y1: ty + 8 }))
        .filter((b) => b.x0 >= L - 30 && b.x1 <= L + pw + Rm);
      tag = opts.find((b) => !hit(b)) || opts[0] || { x0: cur.x + 12, x1: cur.x + 20 + tagW, y0: ty - 8, y1: ty + 8 };
    }
    ctx.fillStyle = C.muted; ctx.textAlign = 'center';
    labels.forEach((l) => {
      if (tag && l.x + 20 > tag.x0 && l.x - 20 < tag.x1 && Math.abs(l.y - (tag.y0 + 8)) < 14) {
        const w = ctx.measureText(l.t).width;
        if (l.x + w / 2 > tag.x0 && l.x - w / 2 < tag.x1) return;
      }
      const w = ctx.measureText(l.t).width;
      roundRect(ctx, l.x - w / 2 - 3, l.y - 7, w + 6, 14, 5); ctx.fillStyle = 'rgba(34,27,42,.88)'; ctx.fill();
      ctx.fillStyle = C.muted; ctx.fillText(l.t, l.x, l.y + 0.5);
    });

    // current marker
    if (d.current) {
      const x = X(((d.current.deg % 360) + 360) % 360), y = Y(d.current.kJ);
      ctx.strokeStyle = 'rgba(255,179,0,.6)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(x + 0.5, T); ctx.lineTo(x + 0.5, T + ph); ctx.stroke(); ctx.setLineDash([]);
      ctx.beginPath(); ctx.arc(x, y, 7, 0, Math.PI * 2);
      ctx.fillStyle = C.mari; ctx.fill(); ctx.lineWidth = 2; ctx.strokeStyle = C.ink; ctx.stroke();
      ctx.font = '700 11px ' + FONT;
      roundRect(ctx, tag.x0, tag.y0, tag.x1 - tag.x0, 16, 8); ctx.fillStyle = 'rgba(15,11,19,.9)'; ctx.fill();
      ctx.fillStyle = C.mari; ctx.textAlign = 'left'; ctx.fillText(tagTxt, tag.x0 + 4, tag.y0 + 8.5);
    }
    ctx.restore();
  }

  /* ---------- geometry helpers ---------- */
  function angleAt(cx, cy, x, y) {
    const a = Math.atan2(x - cx, -(y - cy)) / RAD;
    return (a + 360) % 360;
  }

  // Drag on the annulus 0.5 r – 2.1 r turns the back carbon. onDelta gets the signed
  // clockwise change since the last move, wrapped into (−180, 180].
  function attachDrag(canvas, h) {
    let active = null, last = 0;
    canvas.style.touchAction = 'none';
    const local = (e) => { const b = canvas.getBoundingClientRect(); return [e.clientX - b.left, e.clientY - b.top]; };
    function down(e) {
      if (active != null || (e.pointerType === 'mouse' && e.button !== 0)) return;
      const g = h.getGeometry && h.getGeometry();
      if (!g) return;
      const [x, y] = local(e), dist = Math.hypot(x - g.cx, y - g.cy);
      if (dist < 0.5 * g.r || dist > 2.1 * g.r) return;
      active = e.pointerId; last = angleAt(g.cx, g.cy, x, y);
      try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      e.preventDefault();
      if (h.onStart) h.onStart();
    }
    function move(e) {
      if (e.pointerId !== active) return;
      const g = h.getGeometry(); if (!g) return;
      const [x, y] = local(e), a = angleAt(g.cx, g.cy, x, y);
      let d = a - last; d = ((d + 540) % 360) - 180;
      last = a;
      if (d) h.onDelta(d);
    }
    function up(e) {
      if (e.pointerId !== active) return;
      active = null;
      if (h.onEnd) h.onEnd();
    }
    canvas.addEventListener('pointerdown', down);
    canvas.addEventListener('pointermove', move);
    canvas.addEventListener('pointerup', up);
    canvas.addEventListener('pointercancel', up);
    return function detach() {
      canvas.removeEventListener('pointerdown', down);
      canvas.removeEventListener('pointermove', move);
      canvas.removeEventListener('pointerup', up);
      canvas.removeEventListener('pointercancel', up);
    };
  }

  const api = { draw, drawEnergy, angleAt, attachDrag, chemText, chemWidth, COLORS: C };
  root.NNNewman2D = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

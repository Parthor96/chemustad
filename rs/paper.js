/* rs/paper.js: window.RSPaper, the textbook wedge-and-dash sketch of one stereocenter (canvas 2D).
 * Pure drawing + geometry, node-requirable (tests use a stub context). Contract: rs/SPEC.md section 5.
 *
 * Deviations / additions (documented per the Newman convention):
 *  - geometry(data, ctx?) takes an optional 2D context to measure labels; without one it estimates widths.
 *  - Labels sit just past the bond end (bond end + the label box's half-size along the bond + 6 px)
 *    rather than at exactly 1.25 x bondPx, so wide labels such as CH(CH3)2 never touch their bond.
 *  - drawArrow writes "clockwise"/"counterclockwise" above the drawing (the 90 degree side is always free in
 *    this layout) instead of under the arc's middle, where it would cross the bonds. a.textAt can move it
 *    (degrees, y up), a.text = false hides it.
 *  - PaperData.shiftX (optional) moves the drawing left, so it clears the result card on wide screens.
 *  - PaperData.tag (optional) writes the finished R or S above the drawing.
 *  - Angles everywhere are degrees, counterclockwise from +x with y UP (math convention), as in SLOTS.
 */
(function (root) {
  'use strict';

  const SLOTS = {
    left: { deg: 150, dir: [-0.8165, 0.5774, 0] },
    right: { deg: 30, dir: [0.8165, 0.5774, 0] },
    wedge: { deg: 250, dir: [0, -0.5774, 0.8165] },   // toward the viewer
    dash: { deg: 290, dir: [0, -0.5774, -0.8165] }    // away from the viewer
  };
  const NAMES = ['dash', 'left', 'right', 'wedge'];
  const COL = { night: '#17131d', cream: '#fbf3e2', mari: '#ffb300', teal: '#14b8a6', red: '#e53935', pink: '#ff5c93',
    ink: '#0f0b13', muted: '#c9bfae', dim: '#9d9386' };
  const FONT = 'Inter, -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif';

  const sub = (p, q) => [p[0] - q[0], p[1] - q[1], (p[2] || 0) - (q[2] || 0)];
  function det3(a, b, c) {
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
  }
  const chir = v => Math.sign(det3(sub(v[0], v[3]), sub(v[1], v[3]), sub(v[2], v[3])));

  function perms(a) {
    if (a.length <= 1) return [a.slice()];
    const out = [];
    a.forEach((x, i) => perms(a.slice(0, i).concat(a.slice(i + 1))).forEach(p => out.push([x].concat(p))));
    return out;
  }
  const ALL = perms(NAMES);

  // vecs: bond vectors center -> ligand in rank order (1..4). Returns slot names for ranks 1..4.
  function layout(vecs, pose, variant) {
    const want = chir(vecs);
    const fourOk = pose === 'plane' ? (s => s === 'left' || s === 'right') : (s => s === pose);
    const surv = ALL.filter(p => fourOk(p[3]) && chir(p.map(s => SLOTS[s].dir)) === want);
    surv.sort((x, y) => (x.join(',') < y.join(',') ? -1 : x.join(',') > y.join(',') ? 1 : 0));
    if (!surv.length) return null;
    const k = ((variant | 0) % surv.length + surv.length) % surv.length;
    return surv[k].slice();
  }

  /* ---------------- text helpers (subscript digits), self-contained ---------------- */
  function runs(text) {
    const out = [];
    for (const ch of String(text)) {
      const s = ch >= '0' && ch <= '9';
      const last = out[out.length - 1];
      if (last && last.sub === s) last.t += ch; else out.push({ t: ch, sub: s });
    }
    return out;
  }
  function textWidth(ctx, text, size) {
    const n2 = root.NNNewman2D;
    if (ctx && n2 && n2.chemWidth) return n2.chemWidth(ctx, text, size);
    if (ctx && ctx.measureText) {
      let w = 0;
      runs(text).forEach(r => { ctx.font = '700 ' + (r.sub ? size * 0.72 : size) + 'px ' + FONT; w += ctx.measureText(r.t).width || 0; });
      if (w > 0) return w;
    }
    let w = 0;
    runs(text).forEach(r => { w += r.t.length * (r.sub ? 0.45 : 0.64) * size; });
    return w;
  }
  function drawText(ctx, text, x, y, size, color) {
    const n2 = root.NNNewman2D;
    if (n2 && n2.chemText) return n2.chemText(ctx, text, x, y, size, color);
    const w = textWidth(ctx, text, size);
    let cx = x - w / 2;
    ctx.textBaseline = 'middle'; ctx.textAlign = 'left'; ctx.fillStyle = color;
    runs(text).forEach(r => {
      const s = r.sub ? size * 0.72 : size;
      ctx.font = '700 ' + s + 'px ' + FONT;
      ctx.fillText(r.t, cx, y + (r.sub ? size * 0.23 + 1 : 0));
      cx += ctx.measureText(r.t).width || 0;
    });
    return w;
  }

  /* ---------------- geometry ---------------- */
  const bondPxFor = (w, h) => Math.max(60, Math.min(120, Math.min(w, h) * 0.26));
  const labelSize = w => (w < 420 ? 16 : 18);
  const rad = d => d * Math.PI / 180;
  const sdir = deg => [Math.cos(rad(deg)), -Math.sin(rad(deg))];   // screen direction (y down)

  function geometry(data, ctx) {
    const w = data.w || 300, h = data.h || 300;
    const bondPx = data.bondPx || bondPxFor(w, h);
    const size = labelSize(w);
    const center = [w / 2 - (data.shiftX || 0), h / 2 + bondPx * 0.02];
    const groups = (data.groups || []).map(g => {
      const s = SLOTS[g.slot] || SLOTS.left;
      const d = sdir(s.deg);
      const tw = textWidth(ctx, g.label, size);
      const hw = tw / 2 + 6, hh = size / 2 + 5;
      const edge = Math.min(Math.abs(d[0]) > 1e-6 ? hw / Math.abs(d[0]) : 1e9, Math.abs(d[1]) > 1e-6 ? hh / Math.abs(d[1]) : 1e9);
      const L = bondPx + edge + 4;
      const x = center[0] + d[0] * L, y = center[1] + d[1] * L;
      const outer = Math.min(Math.abs(d[0]) > 1e-6 ? (hw + 13) / Math.abs(d[0]) : 1e9, Math.abs(d[1]) > 1e-6 ? (hh + 13) / Math.abs(d[1]) : 1e9);
      return { x, y, r: Math.max(22, hw + 6), hw, hh, deg: s.deg, end: [center[0] + d[0] * bondPx, center[1] + d[1] * bondPx],
        badge: [x + d[0] * (outer + 1), y + d[1] * (outer + 1)] };
    });
    // long labels on the wedge and the dash would touch: push them apart sideways
    const iw = (data.groups || []).findIndex(g => g.slot === 'wedge'), id = (data.groups || []).findIndex(g => g.slot === 'dash');
    if (iw >= 0 && id >= 0) {
      const a = groups[iw], b = groups[id];
      const over = (a.x + a.hw + 8) - (b.x - b.hw);
      if (over > 0) {
        a.x -= over / 2; a.badge[0] -= over / 2;
        b.x += over / 2; b.badge[0] += over / 2;
      }
    }
    return { center, groups, bondPx, size };
  }

  function hit(data, x, y, ctx) {
    const g = geometry(data, ctx);
    let best = null, bd = Infinity;
    g.groups.forEach((p, i) => {
      const d = Math.hypot(x - p.x, y - p.y);
      if (d <= p.r && d < bd) { bd = d; best = i; }
    });
    return best;
  }

  /* ---------------- drawing ---------------- */
  function drawBadge(ctx, x, y, text, color) {
    ctx.beginPath(); ctx.arc(x, y, 11, 0, Math.PI * 2);
    ctx.fillStyle = color || COL.mari; ctx.fill();
    ctx.lineWidth = 2; ctx.strokeStyle = COL.ink; ctx.stroke();
    ctx.fillStyle = COL.ink; ctx.font = '700 12.5px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText(String(text), x, y + 0.5);
  }

  function drawArrow(ctx, a) {
    const t = a.t == null ? 1 : Math.max(0, Math.min(1, a.t));
    const mod = x => ((x % 360) + 360) % 360;
    const cw = a.dir === 'cw';
    let sweep = cw ? mod(a.from - a.to) : mod(a.to - a.from);
    if (sweep < 1) sweep = 360 - 1e-3;
    const s = sweep * t;
    const color = a.color || COL.mari;
    // canvas angles are clockwise-positive (y down): canvas = -math
    const c0 = rad(-a.from);
    const c1 = cw ? c0 + rad(s) : c0 - rad(s);
    ctx.save();
    ctx.lineCap = 'round';
    ctx.strokeStyle = COL.ink; ctx.lineWidth = 6;
    ctx.beginPath(); ctx.arc(a.cx, a.cy, a.r, c0, c1, !cw); ctx.stroke();
    ctx.strokeStyle = color; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.arc(a.cx, a.cy, a.r, c0, c1, !cw); ctx.stroke();
    if (s > 4) {
      const px = a.cx + a.r * Math.cos(c1), py = a.cy + a.r * Math.sin(c1);
      let tx = -Math.sin(c1), ty = Math.cos(c1);
      if (!cw) { tx = -tx; ty = -ty; }
      const nx = -ty, ny = tx, L = 12, W = 7;
      ctx.beginPath();
      ctx.moveTo(px + tx * L * 0.6, py + ty * L * 0.6);
      ctx.lineTo(px - tx * L * 0.6 + nx * W, py - ty * L * 0.6 + ny * W);
      ctx.lineTo(px - tx * L * 0.6 - nx * W, py - ty * L * 0.6 - ny * W);
      ctx.closePath();
      ctx.fillStyle = color; ctx.fill(); ctx.lineWidth = 1.5; ctx.strokeStyle = COL.ink; ctx.stroke();
    }
    if (a.text !== false && t >= 1) {
      const at = a.textAt == null ? 90 : a.textAt;
      const d = sdir(at), R = a.textR || (a.r + 16);
      ctx.font = '600 13px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      const label = cw ? 'clockwise' : 'counterclockwise';
      const x = a.cx + d[0] * R, y = a.cy + d[1] * R;
      const tw = ctx.measureText ? (ctx.measureText(label).width || 0) : 0;
      ctx.fillStyle = 'rgba(23,19,29,.85)';
      if (ctx.fillRect && tw) ctx.fillRect(x - tw / 2 - 5, y - 9, tw + 10, 18);
      ctx.fillStyle = color; ctx.fillText(label, x, y);
    }
    ctx.restore();
  }

  function draw(ctx, data) {
    const w = data.w, h = data.h;
    const g = geometry(data, ctx);
    const [cx, cy] = g.center;
    ctx.save();
    ctx.fillStyle = COL.night; ctx.fillRect(0, 0, w, h);
    const gap = 14;   // keep bonds clear of the center label
    (data.groups || []).forEach((grp, i) => {
      const p = g.groups[i];
      const d = sdir(p.deg);
      const sx = cx + d[0] * gap, sy = cy + d[1] * gap, ex = p.end[0], ey = p.end[1];
      ctx.globalAlpha = grp.state === 'dim' ? 0.4 : 1;
      if (grp.slot === 'wedge') {
        const nx = -d[1], ny = d[0];
        ctx.beginPath(); ctx.moveTo(sx, sy);
        ctx.lineTo(ex + nx * 4.5, ey + ny * 4.5); ctx.lineTo(ex - nx * 4.5, ey - ny * 4.5); ctx.closePath();
        ctx.fillStyle = COL.cream; ctx.fill();
      } else if (grp.slot === 'dash') {
        const nx = -d[1], ny = d[0];
        ctx.strokeStyle = COL.cream; ctx.lineWidth = 2; ctx.lineCap = 'butt';
        for (let k = 0; k < 7; k++) {
          const f = (k + 0.5) / 7, x = sx + (ex - sx) * f, y = sy + (ey - sy) * f, hw = 1 + 3.5 * f;
          ctx.beginPath(); ctx.moveTo(x + nx * hw, y + ny * hw); ctx.lineTo(x - nx * hw, y - ny * hw); ctx.stroke();
        }
      } else {
        ctx.strokeStyle = COL.cream; ctx.lineWidth = 3; ctx.lineCap = 'round';
        ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(ex, ey); ctx.stroke();
      }
      ctx.globalAlpha = 1;
    });
    // center label
    ctx.globalAlpha = 1;
    // atom number, not a subscript: plain text
    ctx.font = '700 15px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillStyle = COL.cream;
    ctx.fillText(data.centerLabel || 'C', cx, cy);
    // labels
    (data.groups || []).forEach((grp, i) => {
      const p = g.groups[i];
      ctx.globalAlpha = grp.state === 'dim' ? 0.4 : 1;
      if (grp.state === 'ok' || grp.state === 'wrong' || grp.state === 'focus') {
        ctx.beginPath();
        const rx = p.hw + 5, ry = p.hh + 3;
        if (ctx.ellipse) ctx.ellipse(p.x, p.y, rx, ry, 0, 0, Math.PI * 2); else ctx.arc(p.x, p.y, rx, 0, Math.PI * 2);
        ctx.lineWidth = 2.5;
        ctx.strokeStyle = grp.state === 'ok' ? COL.teal : grp.state === 'wrong' ? COL.red : COL.pink;
        ctx.stroke();
      }
      drawText(ctx, grp.label, p.x, p.y, g.size, grp.color || COL.cream);
      if (grp.badge) drawBadge(ctx, p.badge[0], p.badge[1], grp.badge);
      ctx.globalAlpha = 1;
    });
    if (data.arrow) {
      drawArrow(ctx, { cx, cy, r: g.bondPx * 0.55, from: data.arrow.from, via: data.arrow.via, to: data.arrow.to,
        dir: data.arrow.dir, t: data.arrow.t, textAt: 90, textR: g.bondPx * 0.55 + 18, text: data.arrow.text !== false });
    }
    if (data.tag) {
      // R/S result next to the center
      const tx = cx, ty = cy - g.bondPx * 0.55 - 44;
      ctx.font = '400 30px Shrikhand, Georgia, serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillStyle = COL.mari; ctx.fillText(data.tag, tx, Math.max(22, ty));
    }
    if (data.caption) {
      ctx.font = '500 13px ' + FONT; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
      ctx.fillStyle = COL.muted; ctx.fillText(data.caption, cx, h - 14);
    }
    ctx.restore();
  }

  const api = { SLOTS, layout, geometry, draw, hit, drawArrow, drawBadge, bondPxFor };
  root.RSPaper = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

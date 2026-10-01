/* rs/view3d.js: window.RSView3D, the 3D model for the R/S Assigner.
 * three.js (r128 globals) when WebGL works, otherwise a flat 2D-canvas renderer with the same API.
 * Contract: rs/SPEC.md section 6. Written fresh for this page (no code shared with newman/app.js).
 *
 * Frames:
 *  - molecule frame = the coords passed in (Å). local = (coords - centroid) * S, same directions.
 *  - world = local rotated by the model quaternion; the camera looks at a target from a direction.
 *  - screen: CSS px inside the host element, y down.
 *
 * Deviations / additions (documented per the Newman convention):
 *  - create(el, opts) also takes opts.onFrame() (called after every render so the page can draw its
 *    overlay in sync) and returns busy() (true while orient/reset animate).
 *  - orient() takes { front: true } to put `away` toward the camera instead (used by "Put 4 toward me").
 *  - Picking is in screen space for both renderers (nearest projected atom; 22 px on touch, 14 px with a mouse,
 *    or anywhere inside the atom's disc, front-most first) instead of a raycast. One code path serves both.
 *  - setStyle timing (the 900 ms red halo) is handled by the page, which calls setStyle again.
 */
(function (root) {
  'use strict';

  const S = 1.6;
  const DEG = Math.PI / 180;
  const reduced = () => !!(root.matchMedia && root.matchMedia('(prefers-reduced-motion: reduce)').matches);

  const V = {
    sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
    add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
    mul: (a, k) => [a[0] * k, a[1] * k, a[2] * k],
    dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
    cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
    len: a => Math.hypot(a[0], a[1], a[2]),
    unit: a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; },
    lerp: (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t],
    anyPerp: a => { const t = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return V.unit(V.cross(a, t)); }
  };
  // quaternions [x, y, z, w]
  const Q = {
    id: () => [0, 0, 0, 1],
    axis: (ax, ang) => { const s = Math.sin(ang / 2), n = V.unit(ax); return [n[0] * s, n[1] * s, n[2] * s, Math.cos(ang / 2)]; },
    mul: (a, b) => [
      a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
      a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
      a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
      a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]],
    conj: q => [-q[0], -q[1], -q[2], q[3]],
    rot: (q, v) => {
      const x = q[0], y = q[1], z = q[2], w = q[3];
      const tx = 2 * (y * v[2] - z * v[1]), ty = 2 * (z * v[0] - x * v[2]), tz = 2 * (x * v[1] - y * v[0]);
      return [v[0] + w * tx + y * tz - z * ty, v[1] + w * ty + z * tx - x * tz, v[2] + w * tz + x * ty - y * tx];
    },
    // rows of a rotation matrix (world = M * local) -> quaternion
    fromRows: (r0, r1, r2) => {
      const m00 = r0[0], m01 = r0[1], m02 = r0[2], m10 = r1[0], m11 = r1[1], m12 = r1[2], m20 = r2[0], m21 = r2[1], m22 = r2[2];
      const tr = m00 + m11 + m22;
      let x, y, z, w;
      if (tr > 0) { const s = 0.5 / Math.sqrt(tr + 1); w = 0.25 / s; x = (m21 - m12) * s; y = (m02 - m20) * s; z = (m10 - m01) * s; }
      else if (m00 > m11 && m00 > m22) { const s = 2 * Math.sqrt(1 + m00 - m11 - m22); w = (m21 - m12) / s; x = 0.25 * s; y = (m01 + m10) / s; z = (m02 + m20) / s; }
      else if (m11 > m22) { const s = 2 * Math.sqrt(1 + m11 - m00 - m22); w = (m02 - m20) / s; x = (m01 + m10) / s; y = 0.25 * s; z = (m12 + m21) / s; }
      else { const s = 2 * Math.sqrt(1 + m22 - m00 - m11); w = (m10 - m01) / s; x = (m02 + m20) / s; y = (m12 + m21) / s; z = 0.25 * s; }
      const L = Math.hypot(x, y, z, w) || 1;
      return [x / L, y / L, z / L, w / L];
    },
    slerp: (a, b, t) => {
      let bx = b[0], by = b[1], bz = b[2], bw = b[3];
      let c = a[0] * bx + a[1] * by + a[2] * bz + a[3] * bw;
      if (c < 0) { c = -c; bx = -bx; by = -by; bz = -bz; bw = -bw; }
      if (c > 0.9995) {
        const r = [a[0] + (bx - a[0]) * t, a[1] + (by - a[1]) * t, a[2] + (bz - a[2]) * t, a[3] + (bw - a[3]) * t];
        const L = Math.hypot(r[0], r[1], r[2], r[3]); return r.map(v => v / L);
      }
      const th = Math.acos(c), s = Math.sin(th), k0 = Math.sin((1 - t) * th) / s, k1 = Math.sin(t * th) / s;
      return [a[0] * k0 + bx * k1, a[1] * k0 + by * k1, a[2] * k0 + bz * k1, a[3] * k0 + bw * k1];
    }
  };
  const ease = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

  function eigSym3(A) {
    const a = A.map(r => r.slice());
    const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 30; sweep++) {
      if (Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]) < 1e-12) break;
      [[0, 1], [0, 2], [1, 2]].forEach(([p, q]) => {
        if (Math.abs(a[p][q]) < 1e-15) return;
        const th = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]), c = Math.cos(th), s = Math.sin(th);
        for (let k = 0; k < 3; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
        for (let k = 0; k < 3; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
        for (let k = 0; k < 3; k++) { const x = v[k][p], y = v[k][q]; v[k][p] = c * x - s * y; v[k][q] = s * x + c * y; }
      });
    }
    const vals = [a[0][0], a[1][1], a[2][2]];
    return [0, 1, 2].sort((i, j) => vals[j] - vals[i]).map(i => [v[0][i], v[1][i], v[2][i]]);
  }

  const STYLE_COLOR = { center: 0xff5c93, found: 0x14b8a6, wrong: 0xe53935 };
  const EMISSIVE = { ligand: [0xffb300, 0.16], tbA: [0x14b8a6, 0.6], tbB: [0xff5c93, 0.6] };
  const hex = n => '#' + n.toString(16).padStart(6, '0');

  // On the desktop layout the result card floats over the viewer's right side, so the model sits a bit left.
  const viewShift = W => (root.matchMedia && root.matchMedia('(min-width: 861px)').matches ? Math.round(Math.min(110, W * 0.1)) : 0);
  const fitPx = (W, H) => {
    const ox = viewShift(W);
    const right = ox ? W / 2 + ox - 170 : W / 2;
    return 0.8 * Math.max(60, Math.min(H / 2 - 30, W / 2 - ox, right));
  };

  function hasWebGL() {
    if (!root.THREE || !root.THREE.OrbitControls) return false;
    try {
      const c = document.createElement('canvas');
      return !!(root.WebGLRenderingContext && (c.getContext('webgl2') || c.getContext('webgl') || c.getContext('experimental-webgl')));
    } catch (e) { return false; }
  }

  function create(el, opts) {
    opts = opts || {};
    let impl = null;
    if (!opts.forceFlat && hasWebGL()) {
      try { impl = makeGL(el); } catch (e) { console.warn('WebGL failed, using the flat view', e); impl = null; el.innerHTML = ''; }
    }
    if (!impl) impl = makeFlat(el);
    if (!impl) return null;

    // shared state
    const st = { mol: null, coords: null, local: [], centroid: [0, 0, 0], radius: 3, q0: Q.id(), pickable: null, style: {},
      hover: -1, lastKey: '', alive: true };

    function computeLocal() {
      const n = st.coords.length;
      st.centroid = st.coords.reduce((s, c) => V.add(s, V.mul(c, 1 / n)), [0, 0, 0]);
      st.local = st.coords.map(c => V.mul(V.sub(c, st.centroid), S));
    }
    function initialQ() {
      const pts = st.local.filter((p, i) => st.mol.atoms[i].el !== 'H');
      const use = pts.length >= 3 ? pts : st.local;
      const m = [0, 0, 0];
      use.forEach(p => { m[0] += p[0] / use.length; m[1] += p[1] / use.length; m[2] += p[2] / use.length; });
      const cov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      use.forEach(p => { const d = V.sub(p, m); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) cov[i][j] += d[i] * d[j]; });
      const e = eigSym3(cov);
      const x = V.unit(e[0]), y = V.unit(V.sub(e[1], V.mul(x, V.dot(e[1], x)))), z = V.cross(x, y);
      let q = Q.fromRows(x, y, z);
      if (st.mol.rings && st.mol.rings.length) q = Q.mul(Q.axis([1, 0, 0], 0.75), q);
      return q;
    }
    function sceneModel() {
      const chem = root.NNChem, mol = st.mol;
      const nb = mol.atoms.map(() => []);
      mol.bonds.forEach(b => { nb[b.a].push(b.b); nb[b.b].push(b.a); });
      return {
        atoms: mol.atoms.map(a => { const e = (chem && chem.EL[a.el]) || { r: 0.45, color3d: 0x888888 }; return { el: a.el, r: e.r, color: e.color3d }; }),
        bonds: mol.bonds.map(b => ({
          a: b.a, b: b.b, order: b.order || 1,
          toH: mol.atoms[b.a].el === 'H' || mol.atoms[b.b].el === 'H',
          perp: (b.order || 1) > 1 ? pos => {
            const n = nb[b.a].find(j => j !== b.b), m2 = nb[b.b].find(j => j !== b.a);
            if (n != null) return V.sub(pos[n], pos[b.a]);
            if (m2 != null) return V.sub(pos[m2], pos[b.b]);
            return null;
          } : null
        }))
      };
    }

    const view = {
      kind: impl.kind,
      setMolecule(mol, coords) {
        st.mol = mol; st.coords = coords.map(c => c.slice());
        computeLocal();
        st.radius = Math.max(2.2, ...st.local.map((p, i) => V.len(p) + (((root.NNChem && root.NNChem.EL[mol.atoms[i].el]) || { r: 0.45 }).r)));
        st.q0 = initialQ();
        st.style = {}; st.hover = -1; st.pickable = null;
        impl.build(sceneModel(), st.local);
        impl.fit(st.radius, st.q0);
        impl.applyStyle(st.style, st.hover);
      },
      setCoords(coords) {
        st.coords = coords.map(c => c.slice());
        computeLocal();
        impl.setPositions(st.local);
      },
      project(i) {
        if (!st.local[i]) return { x: 0, y: 0, z: 0 };
        return impl.project(st.local[i]);
      },
      projRadius(i) {
        const at = st.mol && st.mol.atoms[i];
        const r = ((root.NNChem && at && root.NNChem.EL[at.el]) || { r: 0.45 }).r * S;
        return impl.pixelSize(st.local[i], r);
      },
      toViewer() { return impl.frame().toViewer; },
      screenFrame() { const f = impl.frame(); return { right: f.right, up: f.up, away: V.mul(f.toViewer, -1) }; },
      orient(o) {
        return new Promise(resolve => {
          if (!st.mol || o.center == null || o.away == null) { resolve(); return; }
          const c = st.coords[o.center];
          let w = V.unit(V.sub(st.coords[o.away], c));            // bond center -> away, molecule frame
          let u;
          if (o.up != null) {
            const d = V.sub(st.coords[o.up], c);
            u = V.sub(d, V.mul(w, V.dot(d, w)));
          }
          if (!u || V.len(u) < 1e-3) u = V.anyPerp(w);
          u = V.unit(u);
          const r2 = o.front ? w : V.mul(w, -1);                   // the molecule vector that maps to world +z (toward camera)
          const r1 = u;
          const r0 = V.cross(r1, r2);
          const q1 = Q.fromRows(r0, r1, r2);
          const focus = st.local[o.center];
          impl.animate(q1, focus, o.ms == null ? 700 : o.ms, resolve, false);
        });
      },
      reset(ms) {
        return new Promise(resolve => { if (!st.mol) { resolve(); return; } impl.animate(st.q0, [0, 0, 0], ms == null ? 700 : ms, resolve, true); });
      },
      setStyle(map) { st.style = map || {}; impl.applyStyle(st.style, st.hover); },
      setPickable(list) { st.pickable = list ? new Set(list) : null; if (st.hover >= 0 && !isPickable(st.hover)) { st.hover = -1; impl.applyStyle(st.style, -1); } },
      resize() { impl.resize(); },
      destroy() { st.alive = false; impl.destroy(); },
      busy: () => impl.busy(),
      get el() { return impl.el; }
    };

    function isPickable(i) {
      if (!st.mol || i < 0) return false;
      if (st.pickable) return st.pickable.has(i);
      return st.mol.atoms[i].el !== 'H';
    }
    function pick(x, y, touch) {
      if (!st.mol) return null;
      const thr = touch ? 22 : 14;
      let best = null, bestZ = -Infinity, near = null, nd = Infinity;
      for (let i = 0; i < st.local.length; i++) {
        if (!isPickable(i)) continue;
        const p = view.project(i), pr = view.projRadius(i);
        const d = Math.hypot(p.x - x, p.y - y);
        if (d <= pr && p.z > bestZ) { bestZ = p.z; best = i; }
        if (d <= Math.max(thr, pr + 4) && d < nd) { nd = d; near = i; }
      }
      return best != null ? best : near;
    }

    // pointers: tap = pointerup with < 6 px movement
    const host = impl.el;
    let down = null;
    const loc = e => { const r = host.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
    host.addEventListener('pointerdown', e => { const p = loc(e); down = { id: e.pointerId, x: p.x, y: p.y, moved: 0, type: e.pointerType }; });
    host.addEventListener('pointermove', e => {
      const p = loc(e);
      if (down && e.pointerId === down.id) { down.moved = Math.max(down.moved, Math.hypot(p.x - down.x, p.y - down.y)); return; }
      if (e.pointerType === 'mouse' && !e.buttons && st.mol && !impl.busy()) {
        const i = pick(p.x, p.y, false);
        const h = i == null ? -1 : i;
        host.style.cursor = h >= 0 ? 'pointer' : '';
        if (h !== st.hover) { st.hover = h; impl.applyStyle(st.style, st.hover); }
      }
    });
    const end = e => {
      if (!down || e.pointerId !== down.id) return;
      const tap = e.type === 'pointerup' && down.moved < 6, type = down.type;
      down = null;
      if (!tap || !st.mol || impl.busy()) return;
      const p = loc(e);
      const i = pick(p.x, p.y, type === 'touch' || type === 'pen');
      if (opts.onTap) opts.onTap(i == null ? null : i);
    };
    host.addEventListener('pointerup', end);
    host.addEventListener('pointercancel', end);
    host.addEventListener('pointerleave', () => { if (st.hover >= 0) { st.hover = -1; host.style.cursor = ''; impl.applyStyle(st.style, -1); } });

    (function loop() {
      if (!st.alive) return;
      requestAnimationFrame(loop);
      impl.render();
      const key = impl.key();
      if (key !== st.lastKey) { st.lastKey = key; if (opts.onChange) opts.onChange(); }
      if (opts.onFrame) opts.onFrame();
    })();

    return view;
  }

  /* ------------------------------------------------------------------ */
  /* three.js renderer                                                    */
  /* ------------------------------------------------------------------ */
  function makeGL(el) {
    const THREE = root.THREE;
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(root.devicePixelRatio || 1, 2));
    el.appendChild(renderer.domElement);
    renderer.domElement.setAttribute('role', 'img');
    renderer.domElement.setAttribute('aria-label', '3D model of the molecule. The atoms and groups are also listed as buttons in the steps panel.');
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 400);
    camera.position.set(0, 2.2, 10);
    const controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true; controls.dampingFactor = 0.08; controls.enablePan = false;
    scene.add(new THREE.AmbientLight(0xffffff, 0.65));
    const d1 = new THREE.DirectionalLight(0xffffff, 0.85); d1.position.set(8, 14, 12); scene.add(d1);
    const d2 = new THREE.DirectionalLight(0xffb300, 0.25); d2.position.set(-10, -4, -6); scene.add(d2);
    const group = new THREE.Group(); scene.add(group);

    let W = 1, H = 1;
    function resize() {
      W = el.clientWidth || 1; H = el.clientHeight || 1;
      renderer.setSize(W, H, false); camera.aspect = W / H;
      const ox = viewShift(W);
      if (ox) camera.setViewOffset(W, H, ox, 0, W, H); else camera.clearViewOffset();
      camera.updateProjectionMatrix();
    }
    const ro = new ResizeObserver(resize); ro.observe(el); resize();

    let model = null, atoms = [], bonds = [], busy = false;
    const halos = {};
    const geo = {};
    const cyl = r => geo['c' + r] || (geo['c' + r] = new THREE.CylinderGeometry(r, r, 1, 18));
    const sph = r => geo['s' + r] || (geo['s' + r] = new THREE.SphereGeometry(r, 32, 32));
    const torus = geo.t || (geo.t = new THREE.TorusGeometry(1, 0.05, 12, 64));
    const camDir0 = V.unit([0, 0.22, 1]);
    let dist0 = 10;

    function clear() {
      while (group.children.length) {
        const o = group.children.pop();
        o.traverse(c => { if (c.material) c.material.dispose(); });
      }
      Object.keys(halos).forEach(k => { scene.remove(halos[k]); halos[k].material.dispose(); delete halos[k]; });
    }
    function build(m, pos) {
      clear(); model = m; atoms = []; bonds = [];
      m.atoms.forEach(a => {
        const mesh = new THREE.Mesh(sph(a.r), new THREE.MeshPhysicalMaterial({ color: a.color, roughness: 0.3, metalness: 0.05, clearcoat: 0.6 }));
        const edge = new THREE.Mesh(sph(+(a.r * 1.06).toFixed(3)), new THREE.MeshBasicMaterial({ color: 0x0f0b13, side: THREE.BackSide }));
        mesh.add(edge);
        group.add(mesh); atoms.push(mesh);
      });
      m.bonds.forEach(b => {
        const g = new THREE.Group();
        const mat = new THREE.MeshStandardMaterial({ color: 0xe8dcc6, roughness: 0.5 });
        const base = b.toH ? 0.09 : 0.12;
        const offs = b.order === 2 ? [-0.14, 0.14] : b.order === 3 ? [-0.2, 0, 0.2] : [0];
        const r = b.order === 1 ? base : b.order === 2 ? 0.085 : 0.075;
        offs.forEach(x => { const c = new THREE.Mesh(cyl(r), mat); c.position.x = x; g.add(c); });
        group.add(g); bonds.push(g);
      });
      setPositions(pos);
    }
    const Yv = new THREE.Vector3(), Xv = new THREE.Vector3(), Zv = new THREE.Vector3(), M = new THREE.Matrix4();
    let curPos = [];
    function setPositions(pos) {
      curPos = pos;
      if (!model) return;
      atoms.forEach((m, i) => m.position.set(pos[i][0], pos[i][1], pos[i][2]));
      model.bonds.forEach((b, k) => {
        const pa = pos[b.a], pb = pos[b.b], d = V.sub(pb, pa), L = V.len(d) || 1e-6;
        const y = V.mul(d, 1 / L);
        let x = b.perp ? b.perp(pos) : null;
        if (!x || V.len(V.sub(x, V.mul(y, V.dot(x, y)))) < 1e-6) x = V.anyPerp(y);
        x = V.unit(V.sub(x, V.mul(y, V.dot(x, y))));
        const z = V.cross(x, y);
        Yv.set(y[0], y[1], y[2]); Xv.set(x[0], x[1], x[2]); Zv.set(z[0], z[1], z[2]);
        M.makeBasis(Xv, Yv, Zv);
        const g = bonds[k];
        g.quaternion.setFromRotationMatrix(M);
        g.position.set((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2, (pa[2] + pb[2]) / 2);
        g.scale.set(1, L, 1);
      });
    }
    function setFade(obj, on) {
      obj.traverse(c => {
        if (!c.material) return;
        if (c.material.opacity === (on ? 0.3 : 1) && c.material.transparent === on) return;
        c.material.transparent = on; c.material.opacity = on ? 0.3 : 1; c.material.depthWrite = !on; c.material.needsUpdate = true;
      });
    }
    function applyStyle(style, hover) {
      if (!model) return;
      atoms.forEach((m, i) => {
        const s = style[i];
        const mat = m.material;
        const e = EMISSIVE[s];
        let col = 0x000000, k = 1;
        if (e) { col = e[0]; k = e[1]; }
        if (i === hover) { if (!e) { col = 0xfff2d0; k = 0.15; } else k += 0.15; }
        mat.emissive.setHex(col); mat.emissiveIntensity = k;
        setFade(m, s === 'dim');
        const hs = STYLE_COLOR[s];
        if (hs != null) {
          let h = halos[i];
          if (!h) { h = new THREE.Mesh(torus, new THREE.MeshBasicMaterial({ color: hs })); scene.add(h); halos[i] = h; }
          h.material.color.setHex(hs);
          const r = model.atoms[i].r * 1.6;
          h.scale.set(r, r, r);
          h.visible = true;
        } else if (halos[i]) halos[i].visible = false;
      });
      bonds.forEach((g, k) => {
        const b = model.bonds[k];
        setFade(g, style[b.a] === 'dim' || style[b.b] === 'dim');
      });
    }
    const tmp = new THREE.Vector3();
    function render() {
      if (controls.enabled) controls.update();
      Object.keys(halos).forEach(k => {
        const h = halos[k]; if (!h.visible || !atoms[k]) return;
        atoms[k].getWorldPosition(tmp); h.position.copy(tmp); h.quaternion.copy(camera.quaternion);
      });
      renderer.render(scene, camera);
    }
    const gq = () => [group.quaternion.x, group.quaternion.y, group.quaternion.z, group.quaternion.w];
    function project(p) {
      const w = Q.rot(gq(), p);
      tmp.set(w[0], w[1], w[2]).project(camera);
      return { x: (tmp.x + 1) / 2 * W, y: (1 - tmp.y) / 2 * H, z: -tmp.z };
    }
    function frame() {
      camera.updateMatrixWorld();
      const inv = Q.conj(gq());
      const qc = [camera.quaternion.x, camera.quaternion.y, camera.quaternion.z, camera.quaternion.w];
      const t = controls.target;
      const v = V.unit([camera.position.x - t.x, camera.position.y - t.y, camera.position.z - t.z]);
      return { toViewer: Q.rot(inv, v), up: Q.rot(inv, Q.rot(qc, [0, 1, 0])), right: Q.rot(inv, Q.rot(qc, [1, 0, 0])) };
    }
    function pixelSize(p, r) {
      if (!p) return 0;
      const a = project(p);
      const f = frame();
      const b = project(V.add(p, V.mul(f.right, r)));
      return Math.hypot(b.x - a.x, b.y - a.y);
    }
    function fit(R, q0) {
      controls.minDistance = Math.max(1.5, 1.3 * R); controls.maxDistance = 5 * R;
      const tanH = Math.tan(22.5 * DEG);
      dist0 = Math.min(controls.maxDistance, Math.max(controls.minDistance, R * (H / 2) / (tanH * fitPx(W, H)) + R * 0.35));
      group.quaternion.set(q0[0], q0[1], q0[2], q0[3]);
      controls.target.set(0, 0, 0);
      camera.position.set(camDir0[0] * dist0, camDir0[1] * dist0, camDir0[2] * dist0);
      camera.up.set(0, 1, 0); camera.lookAt(0, 0, 0);
      controls.update();
    }
    function animate(q1, focusLocal, ms, done, isReset) {
      const qa = gq(), ta = [controls.target.x, controls.target.y, controls.target.z];
      const ca = [camera.position.x, camera.position.y, camera.position.z];
      const d = isReset ? dist0 : Math.min(controls.maxDistance, Math.max(controls.minDistance, camera.position.distanceTo(controls.target)));
      const dir = isReset ? camDir0 : [0, 0, 1];
      const tb = Q.rot(q1, focusLocal), cb = V.add(tb, V.mul(dir, d));
      busy = true; controls.enabled = false;
      const step = t => {
        const q = Q.slerp(qa, q1, t); group.quaternion.set(q[0], q[1], q[2], q[3]);
        const tt = V.lerp(ta, tb, t), cc = V.lerp(ca, cb, t);
        controls.target.set(tt[0], tt[1], tt[2]); camera.position.set(cc[0], cc[1], cc[2]);
        camera.up.set(0, 1, 0); camera.lookAt(controls.target);
      };
      const fin = () => { busy = false; controls.enabled = true; controls.update(); if (done) done(); };
      if (reduced() || ms <= 0) { step(1); fin(); return; }
      const t0 = performance.now();
      (function f() {
        const t = Math.min(1, (performance.now() - t0) / ms);
        step(ease(t));
        if (t < 1) requestAnimationFrame(f); else fin();
      })();
    }
    const r4 = x => Math.round(x * 1e4);
    function key() {
      const q = group.quaternion, p = camera.position, t = controls.target;
      return [q.x, q.y, q.z, q.w, p.x, p.y, p.z, t.x, t.y, t.z, W, H].map(r4).join(',') + '|' + (curPos.length ? curPos[0].map(r4).join(',') : '');
    }
    return {
      kind: 'gl', el: renderer.domElement, build, setPositions, applyStyle, render, project, frame, pixelSize, fit, animate, key,
      resize, busy: () => busy,
      destroy() { ro.disconnect(); clear(); renderer.dispose(); }
    };
  }

  /* ------------------------------------------------------------------ */
  /* flat 2D-canvas renderer (no WebGL)                                   */
  /* ------------------------------------------------------------------ */
  function makeFlat(el) {
    const cv = document.createElement('canvas');
    el.appendChild(cv);
    const ctx = cv.getContext('2d');
    if (!ctx) return null;
    cv.setAttribute('role', 'img');
    cv.setAttribute('aria-label', 'Model of the molecule (flat view). The atoms and groups are also listed as buttons in the steps panel.');
    cv.style.touchAction = 'none'; cv.style.cursor = 'grab';
    let W = 1, H = 1, dpr = 1, q = Q.id(), focus = [0, 0, 0], dist = 10, minD = 3, maxD = 30, dist0 = 10;
    let model = null, pos = [], style = {}, hover = -1, busy = false;
    const tilt = Q.axis([1, 0, 0], 0.22);
    const shade = (n, f) => {
      let r = n >> 16 & 255, g = n >> 8 & 255, b = n & 255;
      const t = f > 0 ? 255 : 0, a = Math.abs(f); r += (t - r) * a; g += (t - g) * a; b += (t - b) * a;
      return 'rgb(' + (r | 0) + ',' + (g | 0) + ',' + (b | 0) + ')';
    };
    function resize() { dpr = Math.min(root.devicePixelRatio || 1, 2); W = el.clientWidth || 1; H = el.clientHeight || 1; cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
    const ro = new ResizeObserver(resize); ro.observe(el); resize();
    function proj(p) {
      const w = Q.rot(q, V.sub(p, focus));
      const f = Math.min(W, H) * 1.1 / Math.max(0.5, dist - w[2]);
      return { x: W / 2 - viewShift(W) + w[0] * f, y: H / 2 - w[1] * f, z: w[2], s: f };
    }
    function render() {
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
      if (!model) return;
      const P = pos.map(proj);
      const list = [];
      model.atoms.forEach((a, i) => list.push({ t: 'a', i, z: P[i].z }));
      model.bonds.forEach((b, k) => list.push({ t: 'b', k, z: (P[b.a].z + P[b.b].z) / 2 - 0.05 }));
      list.sort((u, v) => u.z - v.z);   // far first
      ctx.lineCap = 'round';
      list.forEach(it => {
        if (it.t === 'b') {
          const b = model.bonds[it.k], a = P[b.a], c = P[b.b];
          ctx.globalAlpha = style[b.a] === 'dim' || style[b.b] === 'dim' ? 0.3 : 1;
          const dx = c.x - a.x, dy = c.y - a.y, L = Math.hypot(dx, dy) || 1, nx = -dy / L, ny = dx / L;
          const s = (a.s + c.s) / 2, w = (b.toH ? 0.18 : 0.24) * s;
          const offs = b.order === 2 ? [-0.14, 0.14] : b.order === 3 ? [-0.2, 0, 0.2] : [0];
          const ww = b.order === 1 ? w : w * 0.7;
          offs.forEach(o => {
            const ox = nx * o * s, oy = ny * o * s;
            ctx.strokeStyle = '#0f0b13'; ctx.lineWidth = ww + 3;
            ctx.beginPath(); ctx.moveTo(a.x + ox, a.y + oy); ctx.lineTo(c.x + ox, c.y + oy); ctx.stroke();
            ctx.strokeStyle = '#e8dcc6'; ctx.lineWidth = ww;
            ctx.beginPath(); ctx.moveTo(a.x + ox, a.y + oy); ctx.lineTo(c.x + ox, c.y + oy); ctx.stroke();
          });
        } else {
          const at = model.atoms[it.i], c = P[it.i], r = at.r * c.s, s = style[it.i];
          ctx.globalAlpha = s === 'dim' ? 0.3 : 1;
          const g = ctx.createRadialGradient(c.x - r * 0.35, c.y - r * 0.4, r * 0.1, c.x, c.y, r);
          const lift = it.i === hover ? 0.2 : 0;
          g.addColorStop(0, shade(at.color, 0.55 + lift * 0.5)); g.addColorStop(0.55, shade(at.color, lift)); g.addColorStop(1, shade(at.color, -0.45 + lift));
          ctx.fillStyle = g; ctx.beginPath(); ctx.arc(c.x, c.y, r, 0, 7); ctx.fill();
          const e = EMISSIVE[s];
          ctx.strokeStyle = e ? hex(e[0]) : '#0f0b13'; ctx.lineWidth = e ? (s === 'ligand' ? 2.5 : 3.5) : 2.5; ctx.stroke();
          const hs = STYLE_COLOR[s];
          if (hs != null) { ctx.strokeStyle = hex(hs); ctx.lineWidth = 3; ctx.beginPath(); ctx.arc(c.x, c.y, r * 1.6, 0, 7); ctx.stroke(); }
        }
      });
      ctx.globalAlpha = 1;
    }
    // drag to rotate, pinch / wheel to zoom
    const pts = new Map(); let pinch0 = 0, pd0 = dist, last = null;
    cv.addEventListener('pointerdown', e => {
      if (busy) return;
      try { cv.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      pts.set(e.pointerId, [e.clientX, e.clientY]); last = [e.clientX, e.clientY]; cv.style.cursor = 'grabbing';
      if (pts.size === 2) { const [a, b] = [...pts.values()]; pinch0 = Math.hypot(a[0] - b[0], a[1] - b[1]); pd0 = dist; }
    });
    cv.addEventListener('pointermove', e => {
      if (!pts.has(e.pointerId) || busy) return;
      pts.set(e.pointerId, [e.clientX, e.clientY]);
      if (pts.size === 2) { const [a, b] = [...pts.values()], d = Math.hypot(a[0] - b[0], a[1] - b[1]); if (pinch0) dist = Math.max(minD, Math.min(maxD, pd0 * pinch0 / d)); return; }
      const dx = e.clientX - last[0], dy = e.clientY - last[1]; last = [e.clientX, e.clientY];
      if (Math.hypot(dx, dy) < 0.5) return;
      q = Q.mul(Q.mul(Q.axis([0, 1, 0], dx * 0.01), Q.axis([1, 0, 0], dy * 0.01)), q);
    });
    const up = e => { pts.delete(e.pointerId); if (!pts.size) { last = null; cv.style.cursor = 'grab'; } pinch0 = 0; };
    cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);
    cv.addEventListener('wheel', e => { e.preventDefault(); dist = Math.max(minD, Math.min(maxD, dist * (1 + e.deltaY * 0.001))); }, { passive: false });
    function animate(q1, f1, ms, done, isReset) {
      const qa = q.slice(), fa = focus.slice(), da = dist;
      const qt = isReset ? Q.mul(tilt, q1) : q1, db = isReset ? dist0 : dist;
      busy = true;
      const step = t => { q = Q.slerp(qa, qt, t); focus = V.lerp(fa, f1, t); dist = da + (db - da) * t; };
      const fin = () => { busy = false; if (done) done(); };
      if (reduced() || ms <= 0) { step(1); fin(); return; }
      const t0 = performance.now();
      (function f() {
        const t = Math.min(1, (performance.now() - t0) / ms);
        step(ease(t));
        if (t < 1) requestAnimationFrame(f); else fin();
      })();
    }
    const r4 = x => Math.round(x * 1e4);
    return {
      kind: 'flat', el: cv,
      build(m, p) { model = m; pos = p; },
      setPositions(p) { pos = p; },
      applyStyle(s, h) { style = s || {}; hover = h; },
      render,
      project(p) { const r = proj(p); return { x: r.x, y: r.y, z: r.z }; },
      pixelSize(p, r) { return p ? proj(p).s * r : 0; },
      frame() { const inv = Q.conj(q); return { toViewer: Q.rot(inv, [0, 0, 1]), up: Q.rot(inv, [0, 1, 0]), right: Q.rot(inv, [1, 0, 0]) }; },
      fit(R, q0) { minD = 1.3 * R; maxD = 5 * R; dist0 = dist = Math.max(minD, R * Math.min(W, H) * 1.1 / fitPx(W, H) + R * 0.35); focus = [0, 0, 0]; q = Q.mul(tilt, q0); },
      animate,
      key() { return q.concat(focus, [dist, W, H]).map(r4).join(',') + '|' + (pos.length ? pos[0].map(r4).join(',') : ''); },
      resize, busy: () => busy,
      destroy() { ro.disconnect(); }
    };
  }

  const api = { create, hasWebGL, _Q: Q, _V: V };
  root.RSView3D = api;
})(typeof window !== 'undefined' ? window : globalThis);

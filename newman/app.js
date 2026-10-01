/*
 * newman/app.js  —  page controller for Newman Navigator (newman.html)
 *
 * Wires the input methods (type / draw / presets) to NNChem, builds 3D structures with
 * NNGeom, renders them with three.js (or a flat 2D-canvas renderer when WebGL is not
 * available), and drives bond selection, "look down this bond", the Newman overlay,
 * rotation, the energy plot, and the guided prompts.
 * Contract: newman/SPEC.md section 7. Only global: window.NNApp (test/debug handle).
 *
 * Frames used below:
 *  - "Å frame": st.coords from NNGeom (Ångström).
 *  - "local": the molecule group's frame, local = (coords − centre) × S, S = 1.6 units/Å.
 *    Directions are identical in both frames, so vectors pass straight to NNGeom.
 *  - "world": local rotated by the group quaternion (three.js), or by the flat renderer's
 *    rotation. Screen: CSS px relative to #canvas.
 *
 * Resolved details:
 *  - Picking is done in screen space (nearest projected bond segment, 14 px / 22 px on
 *    touch), for both renderers, instead of a three.js raycast on 0.12-wide cylinders:
 *    thin sticks are hard to hit with a finger, and one code path serves the fallback too.
 *  - The WebGL check uses a throwaway canvas before creating THREE.WebGLRenderer, so a
 *    browser without WebGL never triggers three.js's console.error.
 *  - Overlay hysteresis: it turns on within 12° of the bond axis and off beyond 14°, so it
 *    does not flicker at the edge.
 *  - The side panel always looks from the selected front atom to the back atom with the
 *    front priority group at 12 o'clock (stable while the student rotates); the overlay
 *    matches the 3D view instead, from whichever end the camera is on.
 *  - Dragging the overlay rotates the far half as seen by the viewer (the circle), so the
 *    circle follows the finger from either end.
 *  - The initial orientation puts the molecule's widest spread across the screen
 *    (principal axes of the heavy atoms), then tilts rings so the chair reads as a chair.
 *  - Extra UI beyond the spec: a "Swap ends" button (#swap, spec 7.3 optional) and
 *    clickable examples under the input box.
 *  - Main chain step (CHAIN_SPEC.md 5-6): an optional card between "Choose a molecule" and
 *    "Pick a bond". It never gates the later steps. When the student numbers an equivalent
 *    chain from the other end, setChainOrder() stores it and relabel() redraws every place
 *    that shows C numbers. #chain-msg sits outside #chain-work so "Show me" can report in the
 *    done state. The 3D tag layer is positioned each frame from view.project(). The optional
 *    teal rings on the 3D model during tracing (5.7) are not built.
 */
(function () {
  'use strict';

  const $ = (id) => document.getElementById(id);
  const S = 1.6;                       // scene units per Ångström
  const DEG = Math.PI / 180;
  const ON_DEG = 12, OFF_DEG = 14;     // overlay alignment window
  const reduced = () => !!(window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  const mod360 = (d) => ((d % 360) + 360) % 360;
  const signed = (d) => { const m = mod360(d); return m > 180 ? m - 360 : m; };
  const fmt1 = (x) => (Math.abs(x) < 0.05 ? 0 : x).toFixed(1);
  const coarse = () => !!(window.matchMedia && window.matchMedia('(pointer: coarse)').matches);
  // On the desktop layout the result card and the Newman card float over the right side of
  // the viewer, so the molecule is centred a little left of the viewer's middle.
  const viewShift = (W) => (window.matchMedia && window.matchMedia('(min-width: 861px)').matches ? Math.round(Math.min(110, W * 0.1)) : 0);
  // Largest on-screen radius (px) the molecule may take: ~12% padding inside the part of the
  // viewer left free by the shift (the cards sit on the right).
  const fitPx = (W, H) => {
    const ox = viewShift(W);
    // desktop: keep clear of the cards on the right edge too
    const right = ox ? W / 2 + ox - 180 : W / 2;
    return 0.84 * Math.max(60, Math.min(H / 2, W / 2 - ox, right));
  };
  const esc = (t) => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  // formula / group label -> HTML with subscript digits (CH3 -> CH<sub>3</sub>)
  const subHTML = (t) => esc(t).replace(/([A-Za-z)\]])(\d+)/g, '$1<sub>$2</sub>');

  /* ------------------------------------------------------------------ */
  /* small vector + quaternion helpers (quaternions are [x, y, z, w])     */
  /* ------------------------------------------------------------------ */
  const V3 = {
    sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
    add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
    scale: (a, k) => [a[0] * k, a[1] * k, a[2] * k],
    dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
    cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
    len: (a) => Math.hypot(a[0], a[1], a[2]),
    norm: (a) => { const L = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / L, a[1] / L, a[2] / L]; },
    lerp: (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t],
    perp: (a) => { const t = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return V3.norm(V3.cross(a, t)); }
  };
  const Q = {
    id: () => [0, 0, 0, 1],
    axis: (ax, ang) => { const s = Math.sin(ang / 2), n = V3.norm(ax); return [n[0] * s, n[1] * s, n[2] * s, Math.cos(ang / 2)]; },
    mul: (a, b) => [
      a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
      a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
      a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
      a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]],
    conj: (q) => [-q[0], -q[1], -q[2], q[3]],
    rot: (q, v) => { // q v q*
      const [x, y, z, w] = q, [vx, vy, vz] = v;
      const ix = w * vx + y * vz - z * vy, iy = w * vy + z * vx - x * vz, iz = w * vz + x * vy - y * vx, iw = -x * vx - y * vy - z * vz;
      return [ix * w + iw * -x + iy * -z - iz * -y, iy * w + iw * -y + iz * -x - ix * -z, iz * w + iw * -z + ix * -y - iy * -x];
    },
    // rotation matrix given by its rows (r0, r1, r2) → quaternion
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
      let [bx, by, bz, bw] = b;
      let c = a[0] * bx + a[1] * by + a[2] * bz + a[3] * bw;
      if (c < 0) { c = -c; bx = -bx; by = -by; bz = -bz; bw = -bw; }
      if (c > 0.9995) {
        const r = [a[0] + (bx - a[0]) * t, a[1] + (by - a[1]) * t, a[2] + (bz - a[2]) * t, a[3] + (bw - a[3]) * t];
        const L = Math.hypot(...r); return r.map((v) => v / L);
      }
      const th = Math.acos(c), s = Math.sin(th), k0 = Math.sin((1 - t) * th) / s, k1 = Math.sin(t * th) / s;
      return [a[0] * k0 + bx * k1, a[1] * k0 + by * k1, a[2] * k0 + bz * k1, a[3] * k0 + bw * k1];
    }
  };
  const ease = (t) => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

  // Eigenvectors of a symmetric 3×3 matrix (Jacobi), sorted by eigenvalue, largest first.
  function eigSym3(A) {
    const a = A.map((r) => r.slice());
    let v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 30; sweep++) {
      let off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
      if (off < 1e-12) break;
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        if (Math.abs(a[p][q]) < 1e-15) continue;
        const th = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]);
        const c = Math.cos(th), s = Math.sin(th);
        for (let k = 0; k < 3; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 3; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < 3; k++) { const vkp = v[k][p], vkq = v[k][q]; v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq; }
      }
    }
    const vals = [a[0][0], a[1][1], a[2][2]];
    return [0, 1, 2].sort((i, j) => vals[j] - vals[i]).map((i) => [v[0][i], v[1][i], v[2][i]]);
  }

  /* ------------------------------------------------------------------ */
  /* state                                                                */
  /* ------------------------------------------------------------------ */
  const st = {
    mol: null, coords: null, center: [0, 0, 0], local: [],
    sel: null,              // { bond, front, back, info }
    bondInfos: [], profile: null, labels: null,
    overlayOn: false, overlayNear: -1, overlayFar: -1, overlayGeom: null,
    promptIdx: 0, flat: false, hover: -1, q0: Q.id(), radius: 3,
    dragging: false, dirty: true, lastResultAt: 0, resultPending: false
  };
  let view = null;           // renderer (three.js or flat)
  let sketch = null;         // NNSketch handle (created on first use of Draw)
  const C = () => window.NNChem, G = () => window.NNGeom, N2 = () => window.NNNewman2D;

  /* ------------------------------------------------------------------ */
  /* renderers                                                            */
  /* ------------------------------------------------------------------ */
  function hasWebGL() {
    if (!window.THREE || !THREE.OrbitControls) return false;
    try {
      const c = document.createElement('canvas');
      return !!(window.WebGLRenderingContext && (c.getContext('webgl2') || c.getContext('webgl') || c.getContext('experimental-webgl')));
    } catch (e) { return false; }
  }

  // Shared animation helper: tweens the orientation quaternion and focus point.
  function tween(ms, step, done) {
    if (reduced() || ms <= 0) { step(1); if (done) done(); return; }
    const t0 = performance.now();
    (function f() {
      const t = Math.min(1, (performance.now() - t0) / ms);
      step(ease(t));
      if (t < 1) requestAnimationFrame(f); else if (done) done();
    })();
  }

  /* ----- three.js renderer ----- */
  function makeGL(el) {
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    el.appendChild(renderer.domElement);
    renderer.domElement.setAttribute('role', 'img');
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 400);
    camera.position.set(0, 2.2, 10);
    const controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true; controls.dampingFactor = 0.08; controls.enablePan = false;
    // lights exactly as vsepr.html
    scene.add(new THREE.AmbientLight(0xffffff, 0.65));
    const d1 = new THREE.DirectionalLight(0xffffff, 0.85); d1.position.set(8, 14, 12); scene.add(d1);
    const d2 = new THREE.DirectionalLight(0xffb300, 0.25); d2.position.set(-10, -4, -6); scene.add(d2);
    const group = new THREE.Group(); scene.add(group);
    const halo = new THREE.Mesh(new THREE.TorusGeometry(1, 0.04, 12, 64), new THREE.MeshBasicMaterial({ color: 0xff5c93 }));
    halo.visible = false; scene.add(halo);

    let W = 1, H = 1;
    const resize = () => {
      W = el.clientWidth || 1; H = el.clientHeight || 1;
      renderer.setSize(W, H, false); camera.aspect = W / H;
      const ox = viewShift(W);
      if (ox) camera.setViewOffset(W, H, ox, 0, W, H); else camera.clearViewOffset();
      camera.updateProjectionMatrix();
    };
    new ResizeObserver(resize).observe(el); resize();

    let atoms = [], bonds = [], model = null, style = {}, fadeKey = '';
    const geoCache = {};
    const cyl = (r) => geoCache[r] || (geoCache[r] = new THREE.CylinderGeometry(r, r, 1, 18));
    const sph = (r) => geoCache['s' + r] || (geoCache['s' + r] = new THREE.SphereGeometry(r, 36, 36));
    let camDir0 = V3.norm([0, 0.22, 1]), dist0 = 10;

    function disposeGroup() {
      while (group.children.length) {
        const o = group.children.pop();
        o.traverse((c) => { if (c.material) c.material.dispose(); });
      }
    }
    function build(m) {
      disposeGroup(); model = m; atoms = []; bonds = []; fadeKey = '';
      m.atoms.forEach((a) => {
        const mesh = new THREE.Mesh(sph(a.r), new THREE.MeshPhysicalMaterial({ color: a.color, roughness: 0.3, metalness: 0.05, clearcoat: 0.6 }));
        const edge = new THREE.Mesh(sph(+(a.r * 1.06).toFixed(3)), new THREE.MeshBasicMaterial({ color: 0x0f0b13, side: THREE.BackSide }));
        mesh.add(edge); // ink outline, rickshaw-style
        group.add(mesh); atoms.push(mesh);
      });
      m.bonds.forEach((b, k) => {
        const g = new THREE.Group();
        const mat = new THREE.MeshStandardMaterial({ color: 0xe8dcc6, roughness: 0.5, emissive: 0x000000 });
        const base = b.toH ? 0.09 : 0.12;
        const offs = b.order === 2 ? [-0.14, 0.14] : b.order === 3 ? [-0.2, 0, 0.2] : [0];
        const r = b.order === 1 ? base : b.order === 2 ? 0.085 : 0.075;
        offs.forEach((x) => { const c = new THREE.Mesh(cyl(r), mat); c.position.x = x; g.add(c); });
        g.userData = { bond: k, r };
        group.add(g); bonds.push(g);
      });
      setPositions(m.pos);
    }
    const Y = new THREE.Vector3(), X = new THREE.Vector3(), Z = new THREE.Vector3(), M = new THREE.Matrix4();
    function setPositions(pos) {
      if (!model) return;
      atoms.forEach((m, i) => m.position.set(pos[i][0], pos[i][1], pos[i][2]));
      model.bonds.forEach((b, k) => {
        const pa = pos[b.a], pb = pos[b.b], d = V3.sub(pb, pa), L = V3.len(d) || 1e-6;
        const y = V3.scale(d, 1 / L);
        let x = b.perp ? b.perp(pos) : null;
        if (!x) x = V3.perp(y);
        x = V3.norm(V3.sub(x, V3.scale(y, V3.dot(x, y))));
        const z = V3.cross(x, y);
        Y.set(y[0], y[1], y[2]); X.set(x[0], x[1], x[2]); Z.set(z[0], z[1], z[2]);
        M.makeBasis(X, Y, Z);
        const g = bonds[k];
        g.quaternion.setFromRotationMatrix(M);
        g.position.set((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2, (pa[2] + pb[2]) / 2);
        g.scale.set(1, L, 1);
      });
    }
    function setStyle(s) {
      style = s;
      bonds.forEach((g, k) => {
        const sel = k === s.sel, hov = k === s.hover && !sel;
        const mat = g.children[0].material;
        mat.color.setHex(sel ? 0xffc21a : hov ? 0xfff6e0 : 0xe8dcc6);
        mat.emissive.setHex(sel ? 0xffb300 : hov ? 0x3a2a00 : 0x000000);
        mat.emissiveIntensity = sel ? 0.45 : 1;
        mat.roughness = sel ? 0.35 : 0.5;
        const k2 = sel ? 0.16 / g.userData.r : 1;
        g.children.forEach((c) => { c.scale.x = k2; c.scale.z = k2; });
      });
      const key = s.keep ? Array.from(s.keep.atoms).join(',') + '|' + Array.from(s.keep.bonds).join(',') : '';
      if (key !== fadeKey) {
        fadeKey = key;
        const fade = (obj, on) => obj.traverse((c) => {
          if (!c.material) return;
          c.material.transparent = on; c.material.opacity = on ? 0.3 : 1; c.material.depthWrite = !on; c.material.needsUpdate = true;
        });
        atoms.forEach((m, i) => fade(m, !!s.keep && !s.keep.atoms.has(i)));
        bonds.forEach((g, k) => fade(g, !!s.keep && !s.keep.bonds.has(k)));
      }
      if (s.halo >= 0 && atoms[s.halo]) {
        halo.visible = true;
        const r = model.atoms[s.halo].r * 1.5;
        halo.scale.set(r, r, r);
      } else halo.visible = false;
    }
    const tmp = new THREE.Vector3();
    function frame() {
      if (controls.enabled) controls.update();
      if (halo.visible && atoms[style.halo]) {
        atoms[style.halo].getWorldPosition(tmp); halo.position.copy(tmp); halo.quaternion.copy(camera.quaternion);
      }
      renderer.render(scene, camera);
    }
    const gq = () => [group.quaternion.x, group.quaternion.y, group.quaternion.z, group.quaternion.w];
    function project(p) {
      const w = Q.rot(gq(), p);
      tmp.set(w[0], w[1], w[2]).project(camera);
      return { x: (tmp.x + 1) / 2 * W, y: (1 - tmp.y) / 2 * H, z: tmp.z };
    }
    function basis() {
      camera.updateMatrixWorld();
      const qc = [camera.quaternion.x, camera.quaternion.y, camera.quaternion.z, camera.quaternion.w];
      const inv = Q.conj(gq());
      const t = controls.target;
      const v = V3.norm([camera.position.x - t.x, camera.position.y - t.y, camera.position.z - t.z]);
      return { v: Q.rot(inv, v), up: Q.rot(inv, Q.rot(qc, [0, 1, 0])), right: Q.rot(inv, Q.rot(qc, [1, 0, 0])) };
    }
    function fit(R, q0) {
      controls.minDistance = Math.max(1.5, 1.4 * R); controls.maxDistance = 5 * R;
      const tanH = Math.tan(22.5 * DEG);
      dist0 = Math.min(controls.maxDistance, Math.max(controls.minDistance, R * (H / 2) / (tanH * fitPx(W, H)) + R * 0.35));
      group.quaternion.set(q0[0], q0[1], q0[2], q0[3]);
      controls.target.set(0, 0, 0);
      camera.position.set(camDir0[0] * dist0, camDir0[1] * dist0, camDir0[2] * dist0);
      camera.up.set(0, 1, 0); camera.lookAt(0, 0, 0);
      controls.update();
    }
    // Animate the group to orientation q1 and the camera to look at `focus` (local) from camDir.
    function orient(q1, focusLocal, camDir, dist, instant, done) {
      const qa = gq(), ta = [controls.target.x, controls.target.y, controls.target.z];
      const ca = [camera.position.x, camera.position.y, camera.position.z];
      const d = dist || Math.min(controls.maxDistance, Math.max(controls.minDistance, camera.position.distanceTo(controls.target)));
      const tb = Q.rot(q1, focusLocal), cb = V3.add(tb, V3.scale(camDir, d));
      controls.enabled = false;
      tween(instant ? 0 : 700, (t) => {
        const q = Q.slerp(qa, q1, t); group.quaternion.set(q[0], q[1], q[2], q[3]);
        const tt = V3.lerp(ta, tb, t), cc = V3.lerp(ca, cb, t);
        controls.target.set(tt[0], tt[1], tt[2]); camera.position.set(cc[0], cc[1], cc[2]);
        camera.up.set(0, 1, 0); camera.lookAt(controls.target);
      }, () => { controls.enabled = true; controls.update(); if (done) done(); });
    }
    return {
      kind: 'gl', el: renderer.domElement, build, setPositions, setStyle, frame, project, basis, fit,
      look: (q1, focus, instant, done) => orient(q1, focus, [0, 0, 1], null, instant, done),
      reset: (q0, instant) => orient(q0, [0, 0, 0], camDir0, dist0, instant),
      setControls: (on) => { controls.enabled = on; },
      busy: () => !controls.enabled
    };
  }

  /* ----- flat 2D-canvas renderer (no WebGL) ----- */
  function makeFlat(el) {
    const cv = document.createElement('canvas');
    el.appendChild(cv);
    const ctx = cv.getContext('2d');
    if (!ctx) return null;
    cv.setAttribute('role', 'img');
    cv.style.touchAction = 'none'; cv.style.cursor = 'grab';
    let W = 1, H = 1, dpr = 1, q = Q.id(), focus = [0, 0, 0], dist = 10, minD = 3, maxD = 30;
    let model = null, pos = [], style = {}, busy = false;
    const hex = (n) => '#' + n.toString(16).padStart(6, '0');
    const shade = (n, f) => {
      let r = n >> 16 & 255, g = n >> 8 & 255, b = n & 255;
      const t = f > 0 ? 255 : 0, a = Math.abs(f); r += (t - r) * a; g += (t - g) * a; b += (t - b) * a;
      return 'rgb(' + (r | 0) + ',' + (g | 0) + ',' + (b | 0) + ')';
    };
    function resize() { dpr = Math.min(window.devicePixelRatio || 1, 2); W = el.clientWidth || 1; H = el.clientHeight || 1; cv.width = W * dpr; cv.height = H * dpr; }
    new ResizeObserver(resize).observe(el); resize();
    function project(p) {
      const w = Q.rot(q, V3.sub(p, focus));
      const f = Math.min(W, H) * 1.1 / Math.max(0.5, dist - w[2]);
      return { x: W / 2 - viewShift(W) + w[0] * f, y: H / 2 - w[1] * f, z: -w[2], s: f };
    }
    function frame() {
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
      if (!model) return;
      const P = pos.map(project);
      const list = [];
      model.atoms.forEach((a, i) => list.push({ t: 'a', i, z: P[i].z }));
      model.bonds.forEach((b, k) => list.push({ t: 'b', k, z: (P[b.a].z + P[b.b].z) / 2 + 0.05 }));
      list.sort((u, v) => v.z - u.z); // far first (z is distance-like: bigger = farther)
      ctx.lineCap = 'round';
      const keep = style.keep;
      list.forEach((it) => {
        if (it.t === 'b') {
          const b = model.bonds[it.k], a = P[b.a], c = P[b.b];
          ctx.globalAlpha = keep && !keep.bonds.has(it.k) ? 0.3 : 1;
          const sel = it.k === style.sel, hov = it.k === style.hover;
          const dx = c.x - a.x, dy = c.y - a.y, L = Math.hypot(dx, dy) || 1, nx = -dy / L, ny = dx / L;
          const s = (a.s + c.s) / 2, w = (sel ? 0.32 : b.toH ? 0.18 : 0.24) * s;
          const offs = b.order === 2 ? [-0.14, 0.14] : b.order === 3 ? [-0.2, 0, 0.2] : [0];
          const ww = b.order === 1 ? w : w * 0.7;
          offs.forEach((o) => {
            const ox = nx * o * s, oy = ny * o * s;
            ctx.strokeStyle = '#0f0b13'; ctx.lineWidth = ww + 3;
            ctx.beginPath(); ctx.moveTo(a.x + ox, a.y + oy); ctx.lineTo(c.x + ox, c.y + oy); ctx.stroke();
            ctx.strokeStyle = sel ? '#ffb300' : hov ? '#fff6e0' : '#e8dcc6'; ctx.lineWidth = ww;
            ctx.beginPath(); ctx.moveTo(a.x + ox, a.y + oy); ctx.lineTo(c.x + ox, c.y + oy); ctx.stroke();
          });
        } else {
          const at = model.atoms[it.i], c = P[it.i], r = at.r * c.s;
          ctx.globalAlpha = keep && !keep.atoms.has(it.i) ? 0.3 : 1;
          const g = ctx.createRadialGradient(c.x - r * 0.35, c.y - r * 0.4, r * 0.1, c.x, c.y, r);
          g.addColorStop(0, shade(at.color, 0.55)); g.addColorStop(0.55, hex(at.color)); g.addColorStop(1, shade(at.color, -0.45));
          ctx.fillStyle = g; ctx.beginPath(); ctx.arc(c.x, c.y, r, 0, 7); ctx.fill();
          ctx.strokeStyle = '#0f0b13'; ctx.lineWidth = 2.5; ctx.stroke();
          if (it.i === style.halo) { ctx.strokeStyle = '#ff5c93'; ctx.lineWidth = 3; ctx.beginPath(); ctx.arc(c.x, c.y, r * 1.5, 0, 7); ctx.stroke(); }
        }
      });
      ctx.globalAlpha = 1;
    }
    // drag to rotate, pinch / wheel to zoom
    const pts = new Map(); let pinch0 = 0, dist0 = dist, last = null;
    cv.addEventListener('pointerdown', (e) => {
      if (busy) return;
      try { cv.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      pts.set(e.pointerId, [e.clientX, e.clientY]); last = [e.clientX, e.clientY]; cv.style.cursor = 'grabbing';
      if (pts.size === 2) { const [a, b] = [...pts.values()]; pinch0 = Math.hypot(a[0] - b[0], a[1] - b[1]); dist0 = dist; }
    });
    cv.addEventListener('pointermove', (e) => {
      if (!pts.has(e.pointerId) || busy) return;
      pts.set(e.pointerId, [e.clientX, e.clientY]);
      if (pts.size === 2) { const [a, b] = [...pts.values()], d = Math.hypot(a[0] - b[0], a[1] - b[1]); if (pinch0) dist = Math.max(minD, Math.min(maxD, dist0 * pinch0 / d)); return; }
      const dx = e.clientX - last[0], dy = e.clientY - last[1]; last = [e.clientX, e.clientY];
      q = Q.mul(Q.mul(Q.axis([0, 1, 0], dx * 0.01), Q.axis([1, 0, 0], dy * 0.01)), q);
    });
    const up = (e) => { pts.delete(e.pointerId); if (!pts.size) { last = null; cv.style.cursor = 'grab'; } pinch0 = 0; };
    cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);
    cv.addEventListener('wheel', (e) => { e.preventDefault(); dist = Math.max(minD, Math.min(maxD, dist * (1 + e.deltaY * 0.001))); }, { passive: false });
    function orient(q1, f1, instant, done) {
      const qa = q.slice(), fa = focus.slice();
      busy = true;
      tween(instant ? 0 : 700, (t) => { q = Q.slerp(qa, q1, t); focus = V3.lerp(fa, f1, t); }, () => { busy = false; if (done) done(); });
    }
    const tilt = Q.axis([1, 0, 0], 0.22); // matches the three.js camera's slight elevation
    return {
      kind: 'flat', el: cv,
      build(m) { model = m; pos = m.pos; },
      setPositions(p) { pos = p; },
      setStyle(s) { style = s; },
      frame, project,
      basis() { const inv = Q.conj(q); return { v: Q.rot(inv, [0, 0, 1]), up: Q.rot(inv, [0, 1, 0]), right: Q.rot(inv, [1, 0, 0]) }; },
      fit(R, q0) { minD = 1.4 * R; maxD = 5 * R; dist = Math.max(minD, R * Math.min(W, H) * 1.1 / fitPx(W, H) + R * 0.35); focus = [0, 0, 0]; q = Q.mul(tilt, q0); },
      look(q1, f, instant, done) { orient(q1, f, instant, done); },
      reset(q0, instant) { orient(Q.mul(tilt, q0), [0, 0, 0], instant); },
      setControls() { /* the flat renderer's drag is blocked by the capture listener */ },
      busy: () => busy
    };
  }

  /* ------------------------------------------------------------------ */
  /* molecule → scene                                                     */
  /* ------------------------------------------------------------------ */
  function localPositions() {
    return st.coords.map((c) => V3.scale(V3.sub(c, st.center), S));
  }
  function sceneModel() {
    const mol = st.mol, chem = C();
    const nb = mol.atoms.map(() => []);
    mol.bonds.forEach((b) => { nb[b.a].push(b.b); nb[b.b].push(b.a); });
    return {
      atoms: mol.atoms.map((a) => { const e = chem.EL[a.el] || chem.EL.C; return { el: a.el, r: e.r, color: e.color3d }; }),
      bonds: mol.bonds.map((b) => ({
        a: b.a, b: b.b, order: b.order,
        toH: mol.atoms[b.a].el === 'H' || mol.atoms[b.b].el === 'H',
        // multiple bonds: offset the sticks in the plane of a neighbour
        perp: b.order > 1 ? (pos) => {
          const n = nb[b.a].find((j) => j !== b.b), m = nb[b.b].find((j) => j !== b.a);
          if (n != null) return V3.sub(pos[n], pos[b.a]);
          if (m != null) return V3.sub(pos[m], pos[b.b]);
          return null;
        } : null
      })),
      pos: st.local
    };
  }
  function initialOrientation() {
    const pts = st.local.filter((p, i) => st.mol.atoms[i].el !== 'H');
    const use = pts.length >= 2 ? pts : st.local;
    const cov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const m = [0, 0, 0];
    use.forEach((p) => { m[0] += p[0] / use.length; m[1] += p[1] / use.length; m[2] += p[2] / use.length; });
    use.forEach((p) => { const d = V3.sub(p, m); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) cov[i][j] += d[i] * d[j]; });
    const [e1, e2] = eigSym3(cov);
    const x = V3.norm(e1), y = V3.norm(V3.sub(e2, V3.scale(x, V3.dot(e2, x)))), z = V3.cross(x, y);
    let q = Q.fromRows(x, y, z);
    if (st.mol.rings && st.mol.rings.length) q = Q.mul(Q.axis([1, 0, 0], 0.95), q); // see the chair from the side
    return q;
  }

  /* ------------------------------------------------------------------ */
  /* loading molecules                                                    */
  /* ------------------------------------------------------------------ */
  function say(id, text, kind) { const e = $(id); e.textContent = text || ''; e.className = 'msg' + (kind ? ' ' + kind : ''); }
  function sayHTML(id, html, kind) { const e = $(id); e.innerHTML = html || ''; e.className = 'msg' + (kind ? ' ' + kind : ''); }
  function setPending(on) { const v = document.querySelector('.viewer'); if (v) v.classList.toggle('pending', !!on); }

  function loadInput(text, opts) {
    opts = opts || {};
    const chem = C();
    if (!chem) { say('mol-msg', 'The chemistry module did not load. Reload the page.', 'err'); return false; }
    $('mol-choices').hidden = true;
    setPending(false);
    let r;
    try { r = chem.parse(text); } catch (e) { console.warn('NNChem.parse failed', e); r = { ok: false, error: 'I could not read that. Try a name, a condensed formula, or SMILES.' }; }
    if (r.ok) return loadMol(r.mol, Object.assign({ msgId: 'mol-msg' }, opts));
    if (r.choices) {
      say('mol-msg', r.message || 'That fits more than one molecule. Which one?');
      const box = $('mol-choices');
      box.innerHTML = '';
      r.choices.forEach((ch) => {
        const b = document.createElement('button');
        b.type = 'button'; b.className = 'chip'; b.textContent = ch.label;
        b.addEventListener('click', () => { $('mol-input').value = ch.value; loadInput(ch.value); });
        box.appendChild(b);
      });
      box.hidden = false;
      setPending(true);
      return false;
    }
    say(opts.msgId || 'mol-msg', r.error || 'I could not read that.', 'err');
    return false;
  }

  function loadMol(mol, opts) {
    opts = opts || {};
    const geom = G(), chem = C();
    let built;
    setPending(false);
    try { built = geom.build3D(mol, opts.axial ? { axial: true } : undefined); } catch (e) {
      console.warn('NNGeom.build3D failed', e);
      say(opts.msgId || 'mol-msg', 'I could not build a 3D model of that molecule. Try a smaller or simpler one.', 'err');
      return false;
    }
    st.mol = built.mol; st.coords = built.coords; st.src = mol; st.axial = !!opts.axial;
    // a numbering the student already chose (kept on the source mol) carries over to the rebuild
    if (mol.chainOrder) { try { chem.setChainOrder(st.mol, mol.chainOrder); } catch (e) { console.warn(e); } }
    const n = st.coords.length;
    st.center = st.coords.reduce((s, c) => V3.add(s, V3.scale(c, 1 / n)), [0, 0, 0]);
    st.local = localPositions();
    st.radius = Math.max(2.2, ...st.local.map((p, i) => V3.len(p) + (chem.EL[st.mol.atoms[i].el] || chem.EL.C).r));
    st.sel = null; st.profile = null; st.overlayOn = false; st.hover = -1;
    st.q0 = initialOrientation();
    if (view) { view.build(sceneModel()); view.fit(st.radius, st.q0); }

    try { st.bondInfos = geom.selectableBonds(st.mol); } catch (e) { console.warn(e); st.bondInfos = []; }
    chainLoad(!!opts.keepChain);
    const fml = subHTML(chem.formula(st.mol));
    const shownName = st.mol.name || (st.chain && st.chain.res && st.chain.res.name) || null;
    if (opts.msgId !== false) sayHTML(opts.msgId || 'mol-msg', 'Built ' + (shownName ? esc(shownName) + ' (' + fml + ')' : fml) + '.', 'ok');
    renderBondList();
    renderLegend();
    setCard('c-bond', true);

    // hash for deep links
    if (opts.hash !== false) {
      let h = '';
      if (st.mol.presetId) h = '#' + st.mol.presetId;
      else if (st.mol.source !== 'sketch' && st.mol.input) h = '#smiles=' + encodeURIComponent(st.mol.input);
      try { history.replaceState(null, '', location.pathname + location.search + h); } catch (e) { /* file:// */ }
      lastHash = h;
    }
    if (st.mol.presetId) $('preset').value = st.mol.presetId;

    // choose a bond
    let pick = null;
    const db = st.mol.defaultBond;
    if (opts.keepBond != null && st.bondInfos.some((b) => b.bond === opts.keepBond)) pick = opts.keepBond;
    else if (db) { const k = chem.bondBetween(st.mol, db[0], db[1]); if (st.bondInfos.some((b) => b.bond === k)) pick = k; }
    if (pick == null) { const r = st.bondInfos.find((b) => b.rotatable) || st.bondInfos[0]; if (r) pick = r.bond; }
    st.promptIdx = 0;
    if (pick != null) {
      select(pick);
      if (opts.keepBond == null) startLowest();
    } else {
      say('bond-msg', 'This molecule has no bond with groups on both ends, so there is no Newman projection to draw. Try a longer chain.', 'err');
      ['c-rotate', 'c-energy', 'c-prompts'].forEach((id) => setCard(id, false));
      $('result').hidden = true; $('newman-card').hidden = true; $('rot-strip').hidden = true;
      $('look').disabled = true; $('swap').disabled = true;
      updateStyle();
    }
    return true;
  }

  function setCard(id, open) { const e = $(id); e.classList.toggle('locked', !open); e.inert = !open; }

  // Start an open-chain bond in its lowest-energy anti conformation (else the global
  // minimum), so presets never load in a strained shape.
  function startLowest() {
    if (!st.sel || !st.profile || !st.sel.info.rotatable) return;
    const mins = st.profile.stationary.filter((s) => s.kind === 'min').sort((a, b) => a.kJ - b.kJ);
    const pick = mins.find((s) => s.name === 'anti') || mins[0];
    if (!pick || st.profile.current.kJ <= pick.kJ + 0.05) return;
    st.coords = G().setReferencePhi(st.mol, st.coords, st.sel.front, st.sel.back, pick.deg);
    moved(); refresh(true);
  }

  function renderBondList() {
    const chem = C(), box = $('bond-list');
    box.innerHTML = '';
    const lab = (bi) => chem.bondLabel(st.mol, bi.a, bi.b);
    const key = (bi) => lab(bi).split('–').map((t) => (t[0] === 'C' && /\d/.test(t[1] || '') ? 0 : 100) + (parseInt(t.replace(/\D/g, ''), 10) || 0));
    st.bondInfos.slice().sort((x, y) => { const a = key(x), b = key(y); return a[0] - b[0] || a[1] - b[1]; }).forEach((bi) => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip' + (bi.rotatable ? '' : ' fixed');
      b.dataset.bond = bi.bond;
      b.setAttribute('aria-pressed', 'false');
      b.textContent = chem.bondLabel(st.mol, bi.a, bi.b);
      if (!bi.rotatable) { b.insertAdjacentHTML('afterbegin', LOCK_SVG); b.setAttribute('aria-label', b.textContent + ', locked'); }
      b.addEventListener('click', () => select(bi.bond));
      const tmp = (on) => { if (!st.chain) return; st.chain.tempBond = on ? [bi.a, bi.b] : null; drawPad(); };
      b.addEventListener('mouseenter', () => tmp(true)); b.addEventListener('mouseleave', () => tmp(false));
      b.addEventListener('focus', () => tmp(true)); b.addEventListener('blur', () => tmp(false));
      box.appendChild(b);
    });
    const locked = st.bondInfos.filter((b) => !b.rotatable);
    const ln = $('bond-lock-note');
    ln.hidden = !locked.length;
    if (locked.length) {
      const ring = locked.some((b) => /ring/.test(b.reason || '')), multi = locked.some((b) => !/ring/.test(b.reason || ''));
      ln.innerHTML = LOCK_SVG + ' ' + (ring && multi ? 'Ring bonds and multiple bonds can’t spin.' : ring ? 'Ring bonds can’t spin.' : 'Double and triple bonds can’t spin.');
    }
  }
  const LOCK_SVG = '<svg class="lock" viewBox="0 0 12 14" aria-hidden="true" width="11" height="13"><path d="M3 6V4a3 3 0 0 1 6 0v2" fill="none" stroke="currentColor" stroke-width="1.6"/><rect x="1" y="6" width="10" height="7.5" rx="1.6" fill="currentColor"/></svg>';
  function renderLegend() {
    const chem = C(), seen = [];
    st.mol.atoms.forEach((a) => { if (!seen.includes(a.el)) seen.push(a.el); });
    const order = ['C', 'H', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'B', 'P'];
    seen.sort((a, b) => order.indexOf(a) - order.indexOf(b));
    const tip = coarse() ? 'drag to rotate · pinch to zoom · tap a bond' : 'drag to rotate · scroll to zoom · click a bond';
    $('legend').innerHTML = seen.map((el) => {
      const c = '#' + (chem.EL[el] || chem.EL.C).color3d.toString(16).padStart(6, '0');
      return '<span><i style="background:' + c + ';box-shadow:0 0 0 1px rgba(251,243,226,.45)"></i>' + el + '</span>';
    }).join('') + '<span><i class="lg-bond"></i>selected bond</span><span><i class="lg-front"></i>front atom</span>' +
      '<span class="tip">' + tip + '</span>';
  }

  /* ------------------------------------------------------------------ */
  /* selection                                                            */
  /* ------------------------------------------------------------------ */
  function select(bondIdx, frontOverride) {
    const info = st.bondInfos.find((b) => b.bond === bondIdx);
    if (!info) return false;
    const chem = C(), geom = G();
    let front = Math.min(info.a, info.b), back = Math.max(info.a, info.b);
    if (frontOverride != null && (frontOverride === info.a || frontOverride === info.b)) { front = frontOverride; back = front === info.a ? info.b : info.a; }
    st.sel = { bond: bondIdx, front, back, info };
    if (st.overlayOn) setOverlay(false);
    // cache substituent labels/colours for this bond (the graph does not change while rotating);
    // ring bonds also get an ax / eq tag for each non-ring substituent
    st.labels = {};
    let ae = null;
    if (/ring/.test(info.reason || '')) { try { ae = geom.axialEquatorial(st.mol, st.coords); } catch (e) { ae = null; } }
    [front, back].forEach((side) => {
      chem.neighbors(st.mol, side).forEach((s) => {
        if (s === front || s === back) return;
        const gi = chem.groupInfo(st.mol, s, side);
        st.labels[s] = { label: gi.label, color: (chem.EL[gi.el] || chem.EL.C).color2d, info: gi, tag: ae && ae.get(s) ? ae.get(s) : null };
      });
    });
    try { st.profile = info.rotatable ? geom.torsionProfile(st.mol, st.coords, front, back) : null; } catch (e) { console.warn(e); st.profile = null; }
    let conf0 = null;
    try { conf0 = geom.conformationName(st.mol, st.coords, front, back); } catch (e) { conf0 = null; }

    document.querySelectorAll('#bond-list .chip').forEach((b) => b.setAttribute('aria-pressed', String(+b.dataset.bond === bondIdx)));
    const bl = chem.bondLabel(st.mol, front, back);
    say('bond-msg', info.rotatable ? bl + ' selected. Look down it, then turn the back carbon.' : info.reason, info.rotatable ? 'ok' : '');
    $('look').disabled = false; $('swap').disabled = false;
    ['c-rotate', 'c-energy', 'c-prompts'].forEach((id) => setCard(id, true));

    const rot = info.rotatable;
    const order = (st.mol.bonds[bondIdx] || {}).order || 1;
    const isRing = /ring/.test(info.reason || '');
    $('dihedral').disabled = !rot; $('dihedral-m').disabled = !rot;
    $('rot-controls').hidden = !rot;
    document.querySelectorAll('.snap').forEach((b) => { b.disabled = !rot; });
    $('rot-hint').textContent = rot ? 'Drag the slider, or drag around the Newman circle.'
      : isRing ? 'Locked: ring bond.' : order === 2 ? 'Locked: double bond.' : order === 3 ? 'Locked: triple bond.' : 'Locked: linear end.';
    $('rot-hint').classList.toggle('locked-line', !rot);
    const chair = !!(conf0 && conf0.key === 'chair');
    $('flip').hidden = !chair;
    $('rot-strip').hidden = !rot;
    $('energy-plot').hidden = !rot || !st.profile;
    let note = '';
    if (rot && st.profile) note = energyNote(st.profile);
    else if (rot) note = 'Both ends of this single bond are sp2 (or sp), so conjugation holds it flat (s-trans or s-cis). This simple model does not draw its energy curve.';
    else if (chair) note = chairNote();
    else if (isRing) note = 'Small rings cannot twist their bonds out of line: the C–H bonds on neighbouring carbons stay (nearly) eclipsed, which adds torsional strain to the ring strain.';
    $('energy-note').innerHTML = note;
    $('c-energy').classList.toggle('quiet', !rot);
    $('newman-card').hidden = false;
    $('result').hidden = false;
    st.promptIdx = 0;
    renderPrompt();
    updateStyle();
    refresh(true);
    if (st.chain && st.chain.phase === 'done') { renderMap(); drawPad(); }
    return true;
  }

  // A-values (kJ/mol, approximate): cost of putting a group axial on a cyclohexane chair
  const A_VALUES = { CH3: 7.3, CH2CH3: 7.5, 'CH(CH3)2': 9.2, 'C(CH3)3': 21, OH: 3.9, NH2: 5.9, F: 1.0, Cl: 2.2, Br: 2.0, I: 1.9, SH: 5.0, OCH3: 2.5 };
  function ringSubs() {
    const chem = C(), geom = G(), out = [];
    let ae; try { ae = geom.axialEquatorial(st.mol, st.coords); } catch (e) { return out; }
    ae.forEach((pos, x) => {
      if (st.mol.atoms[x].el === 'H') return;
      const ringAtom = chem.neighbors(st.mol, x).find((k) => (st.mol.rings || []).some((r) => r.includes(k)));
      if (ringAtom == null) return;
      const lab = chem.groupInfo(st.mol, x, ringAtom).label;
      out.push({ atom: x, label: lab, pos, A: A_VALUES[lab] });
    });
    return out;
  }
  function chairNote() {
    const subs = ringSubs();
    if (!subs.length) return 'Every ring carbon has one axial H (straight up or down) and one equatorial H (out around the ring’s edge). Press <b>Flip the chair</b> and they trade places.';
    return subs.map((s) => {
      const g = subHTML(s.label), a = s.A != null ? s.A.toFixed(1) : null;
      if (s.pos === 'eq') return g + ' is equatorial, the lower-energy chair.' + (a ? ' Flip the chair to put it axial: about ' + a + ' kJ/mol higher (its A-value' + (s.label === 'CH3' ? ', two gauche-butane interactions' : '') + ').' : ' Flip the chair to put it axial.');
      return g + ' is axial now, about ' + (a || '?') + ' kJ/mol above the equatorial chair: it bumps the two axial H on the same face (1,3-diaxial strain). Flip back to relieve it.';
    }).join(' ');
  }

  function energyNote(p) {
    if (p.note) return esc(p.note);
    const parts = [];
    (p.terms || []).forEach((t) => {
      const lab = t.label || (t.front + '/' + t.back);
      parts.push(lab + ' eclipsed ' + fmt1(t.ecl));
      if (t.gauche) parts.push(lab + ' gauche ' + fmt1(t.gauche) + (t.oh ? ' (O–H hydrogen bond)' : ''));
    });
    return 'Additive teaching model (approximate kJ/mol): ' + subHTML(parts.join(', ')) + '. Real values shift with solvent.';
  }

  /* ------------------------------------------------------------------ */
  /* rotation                                                             */
  /* ------------------------------------------------------------------ */
  function rotateBy(near, far, delta) {
    if (!st.sel || !st.sel.info.rotatable || !delta) return;
    st.coords = G().rotateAbout(st.mol, st.coords, near, far, delta);
    engage();
    moved();
  }
  function setDihedral(deg) {
    if (!st.sel || !st.sel.info.rotatable) return false;
    st.coords = G().setPriorityDihedral(st.mol, st.coords, st.sel.front, st.sel.back, signed(deg));
    engage();
    moved();
    return true;
  }
  // phone layout: the side Newman card appears once the student has looked or rotated
  function engage() {
    if (st.engaged) return;
    st.engaged = true;
    const v = document.querySelector('.viewer'); if (v) v.classList.add('engaged');
  }
  function moved() {
    st.local = localPositions();
    if (view) view.setPositions(st.local);
    st.dirty = true;
  }

  /* ------------------------------------------------------------------ */
  /* per-frame + UI refresh                                               */
  /* ------------------------------------------------------------------ */
  function newmanFor(near, far, upVec) {
    const chem = C(), geom = G();
    const na = geom.newmanAngles(st.mol, st.coords, near, far, upVec || undefined);
    const pd = geom.dihedralOfPriorityGroups(st.mol, st.coords, st.sel.front, st.sel.back);
    // priority atoms seen from this end: i sits on sel.front, l on sel.back
    const priNear = near === st.sel.front ? pd.i : pd.l, priFar = near === st.sel.front ? pd.l : pd.i;
    const mk = (s) => {
      if (s.lp || s.atom == null) return { angle: s.angle, label: '', color: '#9d9386', lp: true };
      const L = st.labels[s.atom] || { label: chem.atomLabel(st.mol, s.atom), color: '#fbf3e2' };
      return { angle: s.angle, label: L.label, color: L.color, atom: s.atom, tag: L.tag || null };
    };
    const frontSubs = na.frontSubs.map(mk), backSubs = na.backSubs.map(mk);
    const hi = [frontSubs.findIndex((s) => s.atom === priNear), backSubs.findIndex((s) => s.atom === priFar)];
    if (hi[0] >= 0) frontSubs[hi[0]].priority = true;
    if (hi[1] >= 0) backSubs[hi[1]].priority = true;
    return { frontSubs, backSubs, highlight: hi, dihedral: pd.deg, pd };
  }

  function sizeCanvas(cv) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = cv.clientWidth || 1, h = cv.clientHeight || 1;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }

  function panelGeometry() {
    const cv = $('newman2d'), w = cv.clientWidth || 220, h = cv.clientHeight || 220;
    return { cx: w / 2, cy: h / 2, r: Math.max(30, (Math.min(w, h) / 2 - 26) / 1.62) };
  }

  const shortConf = (conf) => (conf && conf.text ? conf.text.split('(')[0].trim() : '');

  // Everything that depends on the current coordinates. `force` skips the 250 ms throttle.
  let lastComputed = null;
  function refresh(force) {
    st.dirty = false;
    if (!st.sel) return;
    const chem = C(), geom = G(), n2 = N2();
    const { front, back, info } = st.sel;
    const nm = newmanFor(front, back, null);
    const pd = nm.pd;
    let conf = null;
    try { conf = geom.conformationName(st.mol, st.coords, front, back); } catch (e) { console.warn(e); }
    const deg = signed(pd.deg);
    let energy = null, es = null;
    if (info.rotatable && st.profile) {
      try { es = geom.energyState(st.mol, st.coords, front, back); energy = es ? es.kJ : null; } catch (e) { console.warn(e); }
    }
    const A = (st.labels[pd.i] || {}).label || 'H', B = (st.labels[pd.l] || {}).label || 'H';
    lastComputed = { pd, conf, energy, A, B, deg };

    // slider(s) + output: signed dihedral, the same number the card and the arc show
    const sv = degInt(deg);
    ['dihedral', 'dihedral-m'].forEach((id) => {
      const sl = $(id);
      if (!(sliderActive && document.activeElement === sl)) { sl.value = sv; if (sliderActive !== id) sliderLast[id] = sv; }
      sl.setAttribute('aria-valuetext', fmtDeg(sv) + ' degrees' + (conf ? ', ' + shortConf(conf).toLowerCase() : ''));
    });
    $('dih-val').textContent = fmtDeg(sv) + '°';
    $('dih-val-m').textContent = fmtDeg(sv) + '°';
    document.querySelectorAll('.snap').forEach((b) => {
      const t = +b.dataset.deg, on = info.rotatable && Math.abs(signed(sv - t)) <= 2;
      b.setAttribute('aria-pressed', String(on));
    });

    // Newman panel
    if (n2) {
      const { ctx, w, h } = sizeCanvas($('newman2d'));
      ctx.clearRect(0, 0, w, h);
      const g = panelGeometry();
      n2.draw(ctx, {
        cx: g.cx, cy: g.cy, r: g.r,
        front: { label: chem.atomLabel(st.mol, front) }, back: { label: chem.atomLabel(st.mol, back) },
        frontSubs: nm.frontSubs, backSubs: nm.backSubs, highlight: nm.highlight,
        dihedral: pd.deg, style: 'panel', rotatable: info.rotatable
      });
      const fl = chem.atomLabel(st.mol, front), bl = chem.atomLabel(st.mol, back);
      $('newman-cap').textContent = 'Looking down ' + fl + ' → ' + bl;
      $('newman2d').setAttribute('aria-label', newmanAria(fl, bl, nm, conf));

      // energy plot: x axis is the reference pair's dihedral (the profile's own angle)
      if (st.profile && info.rotatable && es) {
        const ep = sizeCanvas($('energy-plot'));
        const pri = st.profile.priority || {};
        const la = (st.labels[pri.i] || {}).label || A, lb = (st.labels[pri.l] || {}).label || B;
        n2.drawEnergy(ep.ctx, {
          w: ep.w, h: ep.h, points: st.profile.points, current: { deg: es.deg, kJ: es.kJ },
          stationary: st.profile.stationary,
          xLabel: 'dihedral ' + la + '–' + fl + '–' + bl + '–' + lb, yLabel: 'energy (kJ/mol)'
        });
        $('energy-plot').setAttribute('aria-label', energyAria(es.deg, energy));
      }
    }

    // result card (aria-live): at most every 250 ms while dragging
    const now = performance.now();
    if (force || !st.dragging || now - st.lastResultAt > 250) { writeResult(); st.lastResultAt = now; st.resultPending = false; }
    else st.resultPending = true;
  }
  const fmtDeg = (d) => (d < 0 ? '−' + Math.abs(d) : String(d));
  // one rounding rule for every displayed dihedral: integer in (−180, 180]
  function degInt(d) { const r = Math.round(signed(d)); return r <= -180 ? 180 : r; }
  const RING_KEYS = ['chair', 'boat', 'twist-boat', 'ring', 'ring-eclipsed'];

  function writeResult() {
    if (!st.sel || !lastComputed) return;
    const chem = C(), { pd, conf, energy, A, B } = lastComputed, { front, back, info } = st.sel;
    $('r-conf').textContent = shortConf(conf) || '—';
    $('r-bond').textContent = chem.bondLabel(st.mol, front, back);
    if (conf && conf.ring) {
      const ta = conf.torsionAtoms;
      const lab = ta ? ' (' + ta.map((k) => chem.atomLabel(st.mol, k)).join('–') + ')' : '';
      $('r-dih').textContent = fmtDeg(degInt(conf.dihedral)) + '°' + lab;
    } else {
      $('r-dih').innerHTML = fmtDeg(degInt(pd.deg)) + '° (' + subHTML(A) + ' to ' + subHTML(B) + ')';
    }
    let et;
    if (info.rotatable && energy != null) et = esc(fmt1(energy) + ' kJ/mol');
    else if (info.rotatable) et = 'Not modelled (conjugated)';
    else if ((st.mol.bonds[st.sel.bond] || {}).order > 1) et = 'Fixed (multiple bond)';
    else if (conf && conf.key === 'chair') {
      const subs = ringSubs();
      et = subs.length ? subs.map((x) => subHTML(x.label) + ' ' + (x.pos === 'ax' ? 'axial' + (x.A ? ', +' + x.A.toFixed(1) + ' kJ/mol' : '') : 'equatorial')).join('; ') : 'Chair: no free spin';
    } else if (conf && conf.ring) et = 'Ring bond, no free spin';
    else et = 'Fixed';
    $('r-energy').innerHTML = et;
  }

  function newmanAria(fl, bl, nm, conf) {
    const clock = (a) => { const h = Math.round(a / 30) % 12; return (h === 0 ? 12 : h) + " o'clock"; };
    const list = (subs, withAngle) => subs.map((s, k) => {
      const name = s.lp ? 'lone pair' : s.label;
      if (s.priority) return name + (withAngle ? ' at ' + Math.round(s.angle) + ' degrees' : ' at ' + clock(s.angle));
      return name;
    }).join(', ');
    return 'Newman projection looking down ' + fl + ' to ' + bl + '. Front: ' + list(nm.frontSubs, false) + '. Back: ' + list(nm.backSubs, true) + '. ' + (conf ? shortConf(conf) + '.' : '');
  }
  function energyAria(deg, e) {
    const p = st.profile;
    if (!p) return 'Torsional energy plot';
    const mx = p.stationary.filter((s) => s.kind === 'max').sort((a, b) => b.kJ - a.kJ)[0];
    const mn = p.stationary.filter((s) => s.kind === 'min').sort((a, b) => a.kJ - b.kJ)[0];
    return 'Torsional energy versus dihedral angle. Lowest 0 kJ/mol' + (mn ? ' at ' + mn.name + ' (' + mn.deg + ' degrees)' : '') +
      ', highest ' + fmt1(p.max) + ' kJ/mol' + (mx ? ' at ' + mx.name + ' (' + mx.deg + ' degrees)' : '') +
      '. Now at ' + Math.round(deg) + ' degrees, ' + fmt1(e || 0) + ' kJ/mol.';
  }

  /* ----- style (selection, hover, fade) ----- */
  function updateStyle() {
    if (!view) return;
    let keep = null;
    if (st.overlayOn && st.sel) {
      // Newman view: the drawing takes over; only the two carbons and the sighted bond stay solid
      keep = { atoms: new Set([st.sel.front, st.sel.back]), bonds: new Set([st.sel.bond]) };
    }
    view.setStyle({ sel: st.sel ? st.sel.bond : -1, hover: st.hover, halo: st.sel && !st.overlayOn ? st.sel.front : -1, keep });
  }

  /* ----- overlay ----- */
  function checkOverlay() {
    if (!view || !st.sel) { if (st.overlayOn) setOverlay(false); return; }
    const bs = view.basis();
    const a = V3.norm(V3.sub(st.local[st.sel.back], st.local[st.sel.front]));
    const cf = V3.dot(bs.v, V3.scale(a, -1)); // camera on the front side
    const lim = Math.cos((st.overlayOn ? OFF_DEG : ON_DEG) * DEG);
    let near = -1, far = -1;
    if (cf >= lim) { near = st.sel.front; far = st.sel.back; }
    else if (-cf >= lim) { near = st.sel.back; far = st.sel.front; }
    const on = near >= 0;
    if (on !== st.overlayOn || near !== st.overlayNear) {
      st.overlayNear = near; st.overlayFar = far;
      setOverlay(on);
    }
    if (on) drawOverlay(bs);
  }
  function setOverlay(on) {
    if (on && !st.overlayOn) st.overlayT0 = performance.now();
    st.overlayOn = on;
    $('align-hint').hidden = !on;
    document.querySelector('.viewer').classList.toggle('fading', on);
    if (!on) {
      st.overlayGeom = null;
      const cv = $('newman-overlay'), ctx = cv.getContext('2d');
      ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, cv.width, cv.height);
    }
    updateStyle();
  }
  function drawOverlay(bs) {
    const cv = $('newman-overlay');
    const { ctx, w, h } = sizeCanvas(cv);
    ctx.clearRect(0, 0, w, h);
    const near = st.overlayNear, far = st.overlayFar;
    const pc = view.project(st.local[near]);
    const pr = view.project(V3.add(st.local[near], V3.scale(bs.right, 0.75 * S)));
    const r = Math.max(40, Math.min(140, Math.hypot(pr.x - pc.x, pr.y - pc.y)));
    // overlay canvas and #canvas may not share an origin (phone layout)
    const rc = $('canvas').getBoundingClientRect(), ro = cv.getBoundingClientRect();
    const cx = pc.x + rc.left - ro.left, cy = pc.y + rc.top - ro.top;
    st.overlayGeom = { cx: pc.x, cy: pc.y, r };
    const nm = newmanFor(near, far, bs.up);
    // "the circle appears": a 200 ms fade and scale-in on the first frames
    const t = reduced() ? 1 : Math.min(1, (performance.now() - (st.overlayT0 || 0)) / 200);
    const k = 1 - Math.pow(1 - t, 3);
    N2().draw(ctx, {
      cx, cy, r, front: { label: C().atomLabel(st.mol, near) }, back: { label: C().atomLabel(st.mol, far) },
      frontSubs: nm.frontSubs, backSubs: nm.backSubs, highlight: nm.highlight, dihedral: nm.dihedral,
      style: 'overlay', rotatable: st.sel.info.rotatable, alpha: k, scale: 0.8 + 0.2 * k
    });
  }

  /* ----- chair flip (rebuild with the substituents on the other ring positions) ----- */
  function flipChair() {
    if (!st.sel || !st.src) return;
    const keep = { bond: st.sel.bond, front: st.sel.front, look: st.overlayOn };
    const mol = st.src;
    if (!loadMol(mol, { axial: !st.axial, msgId: false, hash: false, keepBond: keep.bond, keepChain: true })) return;
    if (st.bondInfos.some((b) => b.bond === keep.bond)) select(keep.bond, keep.front);
    if (keep.look) lookDown(true);
    const subs = ringSubs();
    say('bond-msg', subs.length ? 'Chair flipped: ' + subs.map((x) => x.label + ' ' + (x.pos === 'ax' ? 'axial' : 'equatorial')).join(', ') + '.' : 'Chair flipped: every axial H is now equatorial and every equatorial H axial.', 'ok');
  }

  /* ----- look down / reset ----- */
  function lookDown(instant) {
    if (!st.sel || !view) return;
    const f = st.sel.front, b = st.sel.back;
    const a = V3.norm(V3.sub(st.local[b], st.local[f]));
    const pd = G().dihedralOfPriorityGroups(st.mol, st.coords, f, b);
    let u = V3.sub(st.local[pd.i], st.local[f]);
    u = V3.sub(u, V3.scale(a, V3.dot(u, a)));
    u = V3.len(u) < 1e-6 ? V3.perp(a) : V3.norm(u);
    const zl = V3.scale(a, -1), yl = u, xl = V3.cross(yl, zl);
    const q1 = Q.fromRows(xl, yl, zl); // maps local → world with back−front along −z, priority up
    engage();
    view.look(q1, st.local[f], !!instant || reduced(), () => { st.dirty = true; });
    // phones: the viewer sits above the steps, so bring it back into sight to watch the turn
    if (!instant && window.matchMedia && window.matchMedia('(max-width: 860px)').matches) {
      const v = document.querySelector('.viewer');
      if (v && v.getBoundingClientRect().top < -40) v.scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' });
    }
  }

  /* ------------------------------------------------------------------ */
  /* pointer handling on the 3D view: taps pick bonds, overlay drags turn */
  /* ------------------------------------------------------------------ */
  function nearestBond(x, y, thr) {
    let best = null;
    st.mol.bonds.forEach((b, k) => {
      const pa = view.project(st.local[b.a]), pb = view.project(st.local[b.b]);
      const dx = pb.x - pa.x, dy = pb.y - pa.y, L2 = dx * dx + dy * dy || 1;
      let t = ((x - pa.x) * dx + (y - pa.y) * dy) / L2; t = Math.max(0.08, Math.min(0.92, t));
      const d = Math.hypot(x - (pa.x + t * dx), y - (pa.y + t * dy));
      const z = pa.z + (pb.z - pa.z) * t;
      if (d > thr) return;
      if (!best || d < best.d - 3 || (Math.abs(d - best.d) <= 3 && z < best.z)) best = { k, d, z };
    });
    return best;
  }
  function setupViewerPointers() {
    const host = $('canvas');
    let down = null, drag = null;
    const loc = (e) => { const r = host.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
    // Capture phase: a press on the Newman overlay's ring turns the far carbon instead of orbiting.
    host.addEventListener('pointerdown', (e) => {
      if (!st.overlayOn || !st.overlayGeom || !st.sel || !st.sel.info.rotatable) return;
      const p = loc(e), g = st.overlayGeom, d = Math.hypot(p.x - g.cx, p.y - g.cy);
      if (d < 0.5 * g.r || d > 2.1 * g.r) return;
      e.stopPropagation(); e.preventDefault();
      view.setControls(false);
      try { host.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      drag = { id: e.pointerId, last: N2().angleAt(g.cx, g.cy, p.x, p.y) };
      st.dragging = true;
    }, true);
    host.addEventListener('pointerdown', (e) => { if (!drag) down = { id: e.pointerId, ...loc(e), moved: 0, type: e.pointerType }; });
    host.addEventListener('pointermove', (e) => {
      if (drag && e.pointerId === drag.id) {
        const g = st.overlayGeom; if (!g) return;
        const p = loc(e), a = N2().angleAt(g.cx, g.cy, p.x, p.y);
        let d = a - drag.last; d = ((d + 540) % 360) - 180; drag.last = a;
        // positive (clockwise on screen) turns the far atom clockwise as the viewer sees it
        if (d) rotateBy(st.overlayNear, st.overlayFar, d);
        return;
      }
      if (down && e.pointerId === down.id) { const p = loc(e); down.moved = Math.max(down.moved, Math.hypot(p.x - down.x, p.y - down.y)); return; }
      if (e.pointerType === 'mouse' && !e.buttons && st.mol) { // hover feedback
        const p = loc(e), nb = nearestBond(p.x, p.y, 10);
        const k = nb && st.bondInfos.some((b) => b.bond === nb.k) ? nb.k : -1;
        host.style.cursor = k >= 0 ? 'pointer' : '';
        if (k !== st.hover) { st.hover = k; updateStyle(); }
      }
    });
    const end = (e) => {
      if (drag && e.pointerId === drag.id) {
        drag = null; st.dragging = false; view.setControls(true);
        if (st.resultPending) writeResult();
        return;
      }
      if (!down || e.pointerId !== down.id) return;
      const tap = e.type === 'pointerup' && down.moved < 6;
      const type = down.type; down = null;
      if (!tap || !st.mol || view.busy()) return;
      const p = loc(e);
      const nb = nearestBond(p.x, p.y, type === 'touch' ? 22 : 14);
      if (!nb) return;
      if (st.bondInfos.some((b) => b.bond === nb.k)) { if (!st.sel || st.sel.bond !== nb.k) select(nb.k); }
      else say('bond-msg', 'Pick a bond between two atoms that each have other groups attached.', 'err');
    };
    host.addEventListener('pointerup', end);
    host.addEventListener('pointercancel', end);
  }

  /* ------------------------------------------------------------------ */
  /* guided prompts                                                       */
  /* ------------------------------------------------------------------ */
  function pctx() {
    const lc = lastComputed || {};
    const bondLabel = st.sel ? C().bondLabel(st.mol, st.sel.front, st.sel.back) : '';
    return {
      mol: st.mol, coords: st.coords, sel: st.sel, conf: lc.conf, dihedralAbs: lc.pd ? Math.abs(lc.pd.deg) : null,
      energy: lc.energy, profile: st.profile, overlayOn: st.overlayOn,
      groupLabels: { A: lc.A, B: lc.B }, bondLabel, rot: !!(st.sel && st.sel.info.rotatable)
    };
  }
  const pairText = (c) => (c.groupLabels.A === c.groupLabels.B ? 'two ' + c.groupLabels.A + ' groups' : c.groupLabels.A + ' and ' + c.groupLabels.B + ' groups');
  const ethaneLike = (c) => !!(c.profile && c.profile.stationary.some((s) => s.name === 'staggered'));
  const alkylProfile = (c) => !!(c.profile && c.profile.kind !== 'allylic');
  const bothHeavy = (c) => {
    if (!c.sel || !lastComputed) return false;
    const pd = lastComputed.pd, li = st.labels[pd.i], ll = st.labels[pd.l];
    return !!(li && ll && li.info.cls !== 'H' && ll.info.cls !== 'H');
  };
  const dihOf = (i, j, k, l) => G().dihedral(st.coords, i, j, k, l);
  const hetero = (s) => st.labels[s] && ['O', 'N', 'F', 'S'].includes(st.labels[s].info.cls);
  function e2Setup() {
    // leaving group X (Cl/Br/I) on one end, H on the other
    const chem = C(), f = st.sel.front, b = st.sel.back;
    for (const [xEnd, hEnd] of [[f, b], [b, f]]) {
      const X = chem.neighbors(st.mol, xEnd).find((s) => s !== hEnd && ['Cl', 'Br', 'I'].includes(st.mol.atoms[s].el));
      const Hs = chem.neighbors(st.mol, hEnd).filter((s) => s !== xEnd && st.mol.atoms[s].el === 'H');
      if (X != null && Hs.length) return { X, xEnd, hEnd, Hs, xEl: st.mol.atoms[X].el };
    }
    return null;
  }
  function brSetup() {
    const chem = C(), f = st.sel.front, b = st.sel.back;
    const find = (end, other, el) => chem.neighbors(st.mol, end).find((s) => s !== other && st.mol.atoms[s].el === el && (el !== 'C' || true));
    const bf = find(f, b, 'Br'), bb = find(b, f, 'Br'), cf = find(f, b, 'C'), cb = find(b, f, 'C');
    return bf != null && bb != null && cf != null && cb != null ? { bf, bb, cf, cb } : null;
  }

  const PROMPTS = [
    {
      id: 'look', input: 'none', applies: (c) => !!c.sel,
      text: (c) => 'Sight straight down ' + c.bondLabel + '. Rotate the model (or use Look down this bond) until the Newman circle appears.',
      check: (c) => c.overlayOn
        ? { ok: true, msg: 'That is the Newman view. The dot in the middle is the front carbon; the circle is the carbon behind it.' }
        : { ok: false, msg: 'Not lined up yet. Turn the model until you look straight along the bond, or press Look down this bond.' }
    },
    {
      id: 'anti', input: 'number', applies: (c) => c.rot && alkylProfile(c) && bothHeavy(c) && !ethaneLike(c),
      text: (c) => 'Rotate until the ' + pairText(c) + ' are anti. What is the dihedral angle?',
      check: (c, v) => {
        if (!c.conf || c.conf.key !== 'anti') return { ok: false, msg: 'Not anti yet: the ' + pairText(c) + ' are ' + Math.round(c.dihedralAbs) + '° apart. Anti means straight across the circle.' };
        if (Math.abs(v - 180) <= 5) return { ok: true, msg: 'Right: 180°. The two groups sit straight across the circle, as far apart as they can get.' };
        return { ok: false, msg: 'The model is anti, but check the number. Groups straight across the circle are 180° apart.' };
      }
    },
    {
      id: 'gauche', input: 'number', applies: (c) => c.rot && alkylProfile(c) && bothHeavy(c) && !ethaneLike(c),
      text: () => 'Find a gauche conformation. Enter its dihedral.',
      check: (c, v) => {
        if (!c.conf || c.conf.key !== 'gauche') return { ok: false, msg: 'The model is not gauche yet: the ' + pairText(c) + ' are ' + Math.round(c.dihedralAbs) + '° apart. Gauche means neighbours on the circle, 60° apart.' };
        if (Math.abs(v - 60) <= 5 || Math.abs(v - 300) <= 5 || Math.abs(v + 60) <= 5) return { ok: true, msg: 'Yes: ' + Math.round(c.dihedralAbs) + '°. Gauche groups are neighbours on the circle; here that costs ' + fmt1(c.energy || 0) + ' kJ/mol.' };
        return { ok: false, msg: 'The model is gauche, but check the number. Neighbouring positions on the circle are 60° apart.' };
      }
    },
    {
      id: 'max', input: 'none', applies: (c) => c.rot && !!c.profile,
      text: () => 'Rotate to the highest-energy conformation. What is it called?',
      check: (c) => {
        const E = c.energy || 0, p = c.profile;
        if (E < p.max - 0.3) return { ok: false, msg: 'Not the top yet: ' + fmt1(E) + ' kJ/mol here, and the maximum is ' + fmt1(p.max) + '. Watch the marker on the energy plot.' };
        const name = shortConf(c.conf).toLowerCase();
        if (c.conf && c.conf.key === 'totally-eclipsed') {
          const who = c.groupLabels.A === c.groupLabels.B ? 'two ' + c.groupLabels.A + ' groups bump' : c.groupLabels.A + ' and ' + c.groupLabels.B + ' groups bump';
          return { ok: true, msg: 'Yes: totally eclipsed, ' + fmt1(E) + ' kJ/mol. The ' + who + ' into each other.' };
        }
        if (c.profile.kind === 'allylic') {
          return { ok: true, msg: 'Yes: ' + name + ', ' + fmt1(E) + ' kJ/mol. The double bond now sits between two bonds on the sp3 carbon instead of lining up with one. Next to a C=C, that is the high point.' };
        }
        return { ok: true, msg: 'Yes: ' + name + ', ' + fmt1(E) + ' kJ/mol. Every bond on the front carbon lines up with one on the back carbon.' };
      }
    },
    {
      id: 'allyl', input: 'none', applies: (c) => c.rot && !!c.profile && c.profile.kind === 'allylic',
      text: () => 'Next to a double bond the ethane rule flips. Rotate to the lowest energy. What lines up with the double bond?',
      check: (c) => {
        if ((c.energy || 0) > 0.3) return { ok: false, msg: 'Not the lowest yet (' + fmt1(c.energy || 0) + ' kJ/mol). Watch the plot: the minima sit where a bond on the sp3 carbon eclipses the double bond.' };
        const heavy = [c.groupLabels.A, c.groupLabels.B].find((l) => l && l[0] !== '=' && l !== 'H') || 'alkyl';
        const what = c.conf && c.conf.key === 'syn' ? 'the ' + heavy + ' group' : 'a C–H bond';
        return { ok: true, msg: 'Right: ' + what + ' eclipses the double bond here. Around a C=C (or C=O), a bond lined up with the π bond is the favoured shape, the opposite of the ethane rule.' };
      }
    },
    {
      id: 'barrier', input: 'number', applies: (c) => c.rot && ethaneLike(c),
      text: () => 'How much energy (kJ/mol) does it take to go from staggered to eclipsed?',
      check: (c, v) => Math.abs(v - c.profile.max) <= 1
        ? { ok: true, msg: 'Right: about ' + fmt1(c.profile.max) + ' kJ/mol. That is the torsional barrier for this bond.' }
        : { ok: false, msg: 'Not quite. Read the energy at an eclipsed angle (0°, 120°, 240°) and subtract the staggered energy (0).' }
    },
    {
      id: 'e2', input: 'none', applies: (c) => c.rot && !!e2Setup(),
      text: () => { const s = e2Setup(); return 'E2 needs an H and the leaving group anti-periplanar. Rotate until an H on ' + C().atomLabel(st.mol, s.hEnd) + ' is anti to ' + s.xEl + '.'; },
      check: () => {
        const s = e2Setup();
        const best = Math.max(...s.Hs.map((h) => Math.abs(dihOf(s.X, s.xEnd, s.hEnd, h))));
        return best >= 170
          ? { ok: true, msg: 'Yes: that H and the ' + s.xEl + ' are anti-periplanar (' + Math.round(best) + '°), lined up for E2. The C–H and C–' + s.xEl + ' bonds lie in one plane on opposite sides.' }
          : { ok: false, msg: 'The closest H is ' + Math.round(best) + '° from the ' + s.xEl + '. Keep turning until one H sits straight across the circle from it.' };
      }
    },
    {
      id: 'glycol', input: 'none',
      applies: (c) => c.rot && !!c.profile && bothHeavy(c) && hetero(lastComputed.pd.i) && hetero(lastComputed.pd.l) &&
        c.profile.stationary.some((s) => s.kind === 'min' && s.kJ < 0.01 && s.name === 'gauche'),
      text: () => 'In this model, which is lower in energy: anti or gauche? Rotate to the lower one.',
      check: (c) => {
        if (!(c.energy <= 0.3 && c.conf && c.conf.key === 'gauche')) return { ok: false, msg: 'Not the lowest yet (' + fmt1(c.energy || 0) + ' kJ/mol). Compare the energy at anti (180°) and at gauche (60°).' };
        const els = [st.mol.atoms[lastComputed.pd.i].el, st.mol.atoms[lastComputed.pd.l].el];
        const oh = els.includes('O') && (st.labels[lastComputed.pd.i].info.hasH || st.labels[lastComputed.pd.l].info.hasH);
        return { ok: true, msg: oh && !els.includes('F')
          ? 'Right: gauche is lower here. An O–H on one carbon can hydrogen-bond to the O on the other, and that only works when they are close (gauche).'
          : 'Right: gauche is lower. This is the gauche effect: a C–H bond donates electron density into the empty σ* orbital of the C–F (or C–O) bond next to it, which works best when the electronegative groups are gauche.' };
      }
    },
    {
      id: 'ring', input: 'number', applies: (c) => !!(c.sel && c.conf && c.conf.key === 'chair'),
      text: (c) => 'Look down ' + c.bondLabel + '. Every C–C bond in a chair is staggered. What is the C–C–C–C dihedral?',
      check: (c, v) => {
        const t = Math.abs(c.conf.dihedral);
        return Math.abs(Math.abs(v) - t) <= 8
          ? { ok: true, msg: 'About ' + Math.round(t) + '°: gauche. That is why cyclohexane has almost no torsional strain.' }
          : { ok: false, msg: 'Not quite. Read the angle between the two ring bonds in the Newman view (the result card shows it).' };
      }
    },
    {
      id: 'flip', input: 'none', applies: (c) => !!(c.sel && c.conf && c.conf.key === 'chair' && ringSubs().length),
      text: () => { const x = ringSubs()[0]; return 'Press Flip the chair. Is the ' + x.label + ' axial or equatorial now, and which chair is lower in energy?'; },
      check: () => {
        const x = ringSubs()[0];
        return x.pos === 'ax'
          ? { ok: true, msg: 'Right: after the flip the ' + x.label + ' is axial' + (x.A ? ', about ' + x.A.toFixed(1) + ' kJ/mol higher' : '') + '. It bumps the axial H atoms two carbons away (1,3-diaxial strain), so the equatorial chair wins.' }
          : { ok: false, msg: 'The ' + x.label + ' is equatorial now, the lower-energy chair. Press Flip the chair and check again.' };
      }
    },
    {
      id: 'strain', input: 'number', applies: (c) => !!(c.sel && c.conf && c.conf.ring && c.conf.key !== 'chair' && c.conf.size && c.conf.size <= 5),
      text: (c) => 'Look down ' + c.bondLabel + '. Is this ring bond staggered or eclipsed? Enter the ' + (c.conf.size <= 4 ? 'H–C–C–H' : 'C–C–C–C') + ' dihedral (the result card shows it).',
      check: (c, v) => {
        const t = Math.abs(c.conf.dihedral);
        if (Math.abs(Math.abs(v) - t) > 8) return { ok: false, msg: 'Not quite. Read the angle between the two highlighted bonds in the Newman view.' };
        return { ok: true, msg: c.conf.size <= 4
          ? 'Right: about ' + Math.round(t) + '°, (nearly) eclipsed. A small ring cannot twist away from eclipsing, so torsional strain adds to its angle strain.'
          : 'Right: about ' + Math.round(t) + '°. Cyclopentane puckers into an envelope to ease some eclipsing, but its bonds stay far from the 60° of a chair.' };
      }
    },
    {
      id: 'meso', input: 'none', applies: (c) => !!(c.sel && brSetup() && c.mol.presetId && /dibromobutane/.test(c.mol.presetId)),
      text: () => 'Rotate until the two Br are anti. Are the two CH3 groups anti too?',
      check: (c) => {
        const s = brSetup(), br = Math.abs(dihOf(s.bf, st.sel.front, st.sel.back, s.bb)), me = Math.abs(dihOf(s.cf, st.sel.front, st.sel.back, s.cb));
        const name = c.mol.name || 'this isomer';
        if (br < 170) return { ok: false, msg: 'Rotate until the Br atoms are anti first (they are ' + Math.round(br) + '° apart now).' };
        return me >= 170
          ? { ok: true, msg: 'Yes. In ' + name + ', anti Br puts the CH3 groups anti as well (' + Math.round(me) + '°).' }
          : { ok: true, msg: 'No: in ' + name + ' they end up gauche (' + Math.round(me) + '°) when the Br atoms are anti.' };
      }
    }
  ];
  function applicable() { const c = pctx(); return PROMPTS.filter((p) => { try { return p.applies(c); } catch (e) { return false; } }); }
  function renderPrompt() {
    const list = applicable();
    const q = $('prompt-q'), inp = $('prompt-in');
    say('prompt-msg', '');
    if (!list.length) { q.textContent = 'Pick a bond to get a question.'; inp.hidden = true; $('prompt-check').disabled = true; $('prompt-next').disabled = true; return; }
    const p = list[st.promptIdx % list.length];
    q.textContent = p.text(pctx());
    q.dataset.id = p.id;
    inp.hidden = p.input !== 'number'; inp.value = '';
    $('prompt-check').disabled = false; $('prompt-next').disabled = list.length < 2;
  }
  function checkPrompt() {
    if (lastComputed === null) return;
    refresh(true);
    const list = applicable(), p = list.find((x) => x.id === $('prompt-q').dataset.id);
    if (!p) { renderPrompt(); return; }
    let v = null;
    if (p.input === 'number') {
      v = parseFloat($('prompt-in').value);
      if (!isFinite(v)) { say('prompt-msg', 'Type a number of degrees (or kJ/mol) first.', 'err'); $('prompt-in').focus(); return; }
    }
    const r = p.check(pctx(), v);
    say('prompt-msg', r.msg, r.ok ? 'ok' : 'err');
  }

  /* ------------------------------------------------------------------ */
  /* main chain step (CHAIN_SPEC.md sections 5 and 6)                     */
  /* ------------------------------------------------------------------ */
  const SKIP_KEY = 'nn.chainSkip';
  let skipThisLoad = false;   // used when localStorage throws
  const P2 = () => window.NNChain2D;
  const sameArr = (a, b) => !!a && !!b && a.length === b.length && a.every((x, k) => x === b[k]);
  const STEREO_RE = /^(\((?:\d?[RSEZ],?)+\)-|meso-|cis-|trans-)/i;
  function skipRemembered() {
    try { if (window.localStorage && localStorage.getItem(SKIP_KEY) === '1') return true; } catch (e) { /* storage blocked */ }
    return skipThisLoad;
  }
  function rememberSkip(on) {
    try { if (on) localStorage.setItem(SKIP_KEY, '1'); else localStorage.removeItem(SKIP_KEY); skipThisLoad = false; }
    catch (e) { skipThisLoad = !!on; }
  }
  function newChainState() {
    return { phase: 'na', res: null, path: [], fails: 0, shorterFails: 0, ringC1: null, order: null, tags: false,
      pad: st.chain ? st.chain.pad : null, layout: null, flash: null, numbersDim: null, highlight: [], sticky: null,
      tempBond: null, anim: 0, showing: false };
  }
  st.chain = newChainState();
  const ch = () => st.chain;
  const pcgLabel = (res) => (res && res.pcg ? res.pcg.label : 'OH');

  // Called from loadMol after st.mol is set. keep = true keeps the phase (chair flip rebuilds).
  function chainLoad(keep) {
    const chem = C();
    const old = st.chain;
    if (keep && old.res) {
      try { old.res = chem.mainChain(st.mol); } catch (e) { console.warn(e); }
      if (old.order) { chem.setChainOrder(st.mol, old.order); }
      if (old.pad && old.layout) old.pad.setMol(st.mol, old.layout);
      renderChain();
      return;
    }
    const s = newChainState();
    st.chain = s;
    s.anim++;
    try {
      s.res = chem.mainChain(st.mol);
      const k = s.res.kind;
      if (!P2() || k === 'none' || !s.res.chain.length) s.phase = 'na';
      else if (k === 'single') { s.phase = 'done'; s.tags = true; }
      else if (skipRemembered()) s.phase = 'skipped';
      else s.phase = 'intro';
    } catch (e) {
      console.warn('main chain step failed', e);
      s.res = null; s.phase = 'na';
    }
    setCard('c-chain', true);
    renderChain();
  }
  function ensurePad() {
    const s = ch();
    if (!P2()) return null;
    $('chain-wrap').hidden = false;
    if (!s.pad) {
      s.pad = P2().attach($('chain-pad'), {
        onTap: (i, info) => chainTap(i, info),
        onEndTap: (i) => chainEnd(i),
        onUndo: () => chainUndo()
      });
    }
    if (!s.layout) {
      try { s.layout = P2().layout(st.mol, { src: st.src }); } catch (e) { console.warn('chain layout failed', e); s.layout = null; return null; }
      s.pad.setMol(st.mol, s.layout);
    }
    return s.pad;
  }
  function mainOrder() {
    const s = ch();
    if (!s.res) return [];
    return (st.mol.chainOrder && st.mol.chainOrder.length) ? st.mol.chainOrder : s.res.chain;
  }
  function chainCandidates() {
    const s = ch(), chem = C();
    if (s.phase === 'trace') {
      if (!s.path.length) {
        // before the first tap: ring the end carbons (one carbon neighbor), where a main chain starts
        const inRing = new Set([].concat.apply([], st.mol.rings || []));
        return st.mol.atoms.map((a, i) => i).filter((i) => st.mol.atoms[i].el === 'C' && !inRing.has(i) &&
          chem.neighbors(st.mol, i).filter((j) => st.mol.atoms[j].el === 'C').length <= 1);
      }
      const last = s.path[s.path.length - 1];
      return chem.neighbors(st.mol, last).filter((j) => st.mol.atoms[j].el === 'C' && s.path.indexOf(j) < 0);
    }
    if (s.phase === 'ring' && s.ringC1 != null) {
      const ring = s.res.chain;
      return chem.neighbors(st.mol, s.ringC1).filter((j) => ring.indexOf(j) >= 0);
    }
    return [];
  }
  function padState() {
    const s = ch(), chem = C(), out = { path: [], candidates: [], numbers: {}, numbersDim: {}, branchLabels: {}, ends: [], flash: s.flash, highlight: s.highlight || [], selBond: null, tracing: false };
    if (s.phase === 'trace' || s.phase === 'ring' || s.phase === 'direction') {
      out.path = s.phase === 'ring' ? (s.ringC1 != null ? [s.ringC1] : []) : s.path.slice();
      out.candidates = chainCandidates();
      out.tracing = true;
      if (s.phase === 'direction') out.ends = [s.path[0], s.path[s.path.length - 1]];
      if (s.numbersDim) out.numbersDim = s.numbersDim;
      if (s.showing) { out.path = s.path.slice(); out.candidates = []; out.ends = []; out.numbers = s.showNumbers || {}; }
    } else if (s.phase === 'done') {
      const order = mainOrder();
      order.forEach((a, k) => { out.numbers[a] = k + 1; });
      st.mol.atoms.forEach((a, i) => { if (a.el === 'C' && order.indexOf(i) < 0) out.branchLabels[i] = chem.atomLabel(st.mol, i); });
      const sb = s.tempBond || (st.sel ? [st.sel.front, st.sel.back] : null);
      if (sb && st.mol.atoms[sb[0]].el !== 'H' && st.mol.atoms[sb[1]].el !== 'H') out.selBond = sb;
    }
    return out;
  }
  function drawPad() { const s = ch(); if (s.pad && s.layout && !$('chain-wrap').hidden) s.pad.setState(padState()); }
  function flash(atoms, kind) { ch().flash = atoms && atoms.length ? { atoms: atoms.slice(), kind } : null; }
  function chainSay(text, kind) { say('chain-msg', text, kind); }

  function chainQuestion() {
    const s = ch(), res = s.res;
    if (s.phase === 'trace') {
      let q = 'Tap the carbons of the longest chain, one after another, from one end to the other. You can also drag through them.';
      if (res.pcg && (s.fails > 0 || /^(acid|aldehyde|nitrile|amide)$/.test(res.pcg.cls))) q += ' The main chain has to include the carbon with the ' + res.pcg.label + '.';
      if (res.unsat.length && s.fails > 0) q += ' It also has to include both carbons of the ' + (res.unsat[0].order === 3 ? 'C≡C' : 'C=C') + '.';
      return q;
    }
    if (s.phase === 'ring') return s.ringC1 == null
      ? 'In a ring, the ring is the main chain. Tap the carbon that should be C1, then the one that should be C2.'
      : 'Now tap the carbon next to it that should be C2.';
    if (s.phase === 'direction') return 'Now number it. Tap the end that should be C1.';
    return '';
  }
  function renderChain() {
    const s = ch(), show = (id, on) => { $(id).hidden = !on; };
    const ph = s.phase;
    show('chain-intro', ph === 'intro');
    show('chain-work', ph === 'trace' || ph === 'ring' || ph === 'direction');
    show('chain-done', ph === 'done');
    show('chain-skipped', ph === 'skipped');
    show('chain-na', ph === 'na');
    const padOn = ph === 'trace' || ph === 'ring' || ph === 'direction' || ph === 'done';
    $('chain-wrap').hidden = !padOn;
    if (padOn) ensurePad();
    if (ph === 'na') $('chain-na-text').textContent = s.res && s.res.reasonText && s.res.kind === 'none' ? s.res.reasonText
      : 'This tool numbers single chains and single carbon rings. This molecule has more than that, so I number it the simple way and skip the naming step.';
    if (ph !== 'trace' && ph !== 'ring' && ph !== 'direction' && ph !== 'done') chainSay('');
    $('chain-q').textContent = chainQuestion();
    $('chain-check').hidden = ph !== 'trace';
    $('chain-reset').hidden = ph === 'direction';
    $('chain-undo').disabled = ph === 'trace' ? !s.path.length : ph === 'ring' ? s.ringC1 == null : false;
    $('chain-show').disabled = !!s.showing;
    $('chain-check').disabled = !!s.showing;
    const n = ph === 'trace' ? s.path.length : 0;
    $('chain-count').textContent = n ? n + (n === 1 ? ' carbon' : ' carbons') : '';
    if (ph === 'done') renderDone();
    $('chain-tags').setAttribute('aria-pressed', String(!!s.tags));
    $('chain-tags2').setAttribute('aria-pressed', String(!!s.tags));
    const nm = st.mol && s.res && s.res.name ? s.res.name : (st.mol && st.mol.name) || 'this molecule';
    $('chain-pad').setAttribute('aria-label', 'Skeletal structure of ' + nm + (ph === 'done' ? ', main chain numbered' : '') + '. Arrow keys move between atoms, Enter picks one, Backspace undoes.');
    drawPad();
    buildTags();
  }

  /* ----- done state: name breakdown ----- */
  function typedName() {
    const m = st.src || st.mol;
    if (!m) return null;
    const norm = (t) => String(t || '').toLowerCase().replace(/[\s_]/g, '');
    if (m.source === 'name' && m.input) {
      if (m.presetId && norm(m.input).replace(/-/g, '') === norm(m.presetId).replace(/-/g, '')) return m.name ? { text: m.name, typed: false } : null;
      return { text: m.input.trim(), typed: true };
    }
    return m.name ? { text: m.name, typed: false } : null;
  }
  function renderDone() {
    const s = ch(), res = s.res, chem = C();
    const box = $('chain-name'), list = $('chain-lines');
    box.innerHTML = ''; list.innerHTML = '';
    box.classList.toggle('unnamed', !res.supported);
    if (!res.supported || !res.parts) {
      box.textContent = st.mol.name || '';
      $('chain-alt').textContent = res.reasonText || '';
    } else {
      res.parts.forEach((p, k) => {
        if (p.role === 'punct' || p.line == null) {
          const sp = document.createElement('span'); sp.className = 'part ' + p.role; sp.textContent = p.text; box.appendChild(sp); return;
        }
        const b = document.createElement('button');
        b.type = 'button'; b.className = 'part ' + p.role; b.textContent = p.text;
        b.dataset.line = p.line; b.dataset.k = k;
        b.setAttribute('aria-describedby', 'chain-line-' + p.line);
        b.addEventListener('mouseenter', () => hiLine(p.line, false));
        b.addEventListener('mouseleave', () => hiLine(null, false));
        b.addEventListener('focus', () => hiLine(p.line, false));
        b.addEventListener('blur', () => hiLine(null, false));
        b.addEventListener('click', () => hiLine(p.line, true));
        box.appendChild(b);
      });
      (res.lines || []).forEach((l, k) => {
        const li = document.createElement('li');
        li.id = 'chain-line-' + k; li.tabIndex = 0;
        li.innerHTML = esc(l.text).replace(/([H)])(\d+)/g, '$1<sub>$2</sub>'); // CH3, NH2, C(CH3)3; never C2
        li.addEventListener('mouseenter', () => hiLine(k, false));
        li.addEventListener('mouseleave', () => hiLine(null, false));
        li.addEventListener('focus', () => hiLine(k, false));
        li.addEventListener('blur', () => hiLine(null, false));
        li.addEventListener('click', () => hiLine(k, true));
        list.appendChild(li);
      });
      // typed / common name, or the 2013 form
      const stereo = res.stereo || '';
      const tn = typedName();
      const norm = (t) => String(t || '').toLowerCase().replace(/\s+/g, '').replace(STEREO_RE, '');
      let alt = '';
      if (tn && norm(tn.text) !== norm(res.name) && norm(tn.text) !== norm(res.name2013)) {
        alt = (tn.typed ? 'You typed ' + tn.text + '. ' : 'Common name: ' + tn.text + '. ') + 'The IUPAC name is ' + stereo + res.name + '.';
      } else if (tn && norm(tn.text) === norm(res.name2013) && res.name2013 !== res.name) {
        alt = 'Also written: ' + stereo + res.name + '.';
      } else if (res.name2013 && res.name2013 !== res.name) alt = 'Also written: ' + stereo + res.name2013 + '.';
      if (res.unsat.length && res.maxLen > res.n && res.kind === 'chain') {
        alt += (alt ? ' ' : '') + 'Some newer textbooks pick the longest chain even when it misses the C=C; this tool keeps the C=C in the chain.';
      }
      $('chain-alt').textContent = alt;
    }
    $('chain-alt').hidden = !$('chain-alt').textContent;
    renderMap();
  }
  // connect the numbers to the bond list (#chain-map)
  function renderMap() {
    const s = ch(), chem = C();
    if (!s.res || !st.mol) return;
    let map = '';
    if (st.sel) {
      const bl = chem.bondLabel(st.mol, st.sel.front, st.sel.back), parts = bl.split('–');
      map = bl + ' in the bond list is the bond between ' + parts[0] + ' and ' + parts[1] + '.';
    } else map = 'The bond list uses these numbers.';
    const order = mainOrder();
    const branch = st.mol.atoms.map((a, i) => i).filter((i) => st.mol.atoms[i].el === 'C' && order.indexOf(i) < 0)
      .map((i) => chem.atomLabel(st.mol, i)).sort((a, b) => parseInt(a.slice(1), 10) - parseInt(b.slice(1), 10));
    if (branch.length === 1) map += ' The branch carbon gets ' + branch[0] + ' after the chain.';
    else if (branch.length > 1) map += ' The branch carbons get ' + (branch.length > 2 ? branch.slice(0, -1).join(', ') + ' and ' + branch[branch.length - 1] : branch.join(' and ')) + ' after the chain.';
    if (s.tags) map += ' The 3D model uses the same numbers.';
    $('chain-map').textContent = map;
  }
  function hiLine(k, sticky) {
    const s = ch(), res = s.res;
    if (sticky) s.sticky = s.sticky === k ? null : k;
    const line = k != null ? k : s.sticky;
    s.highlight = line != null && res && res.lines && res.lines[line] ? res.lines[line].atoms.filter((a) => st.mol.atoms[a] && st.mol.atoms[a].el !== 'H') : [];
    document.querySelectorAll('#chain-name .part').forEach((b) => b.classList.toggle('on', line != null && b.dataset.line === String(line)));
    document.querySelectorAll('#chain-lines li').forEach((li, j) => li.classList.toggle('on', j === line));
    drawPad();
  }

  /* ----- actions ----- */
  function chainStart() {
    const s = ch();
    if (!s.res || s.phase === 'na') return false;
    s.anim++; s.showing = false;
    s.path = []; s.ringC1 = null; s.flash = null; s.numbersDim = null; s.fails = 0; s.shorterFails = 0;
    if (s.res.kind === 'single') s.phase = 'done';
    else s.phase = s.res.kind === 'chain' ? 'trace' : 'ring';
    chainSay('');
    renderChain();
    return true;
  }
  function chainSkip() {
    const s = ch();
    s.anim++; s.showing = false;
    rememberSkip(true);
    s.phase = 'skipped'; s.tags = false;
    renderChain();
  }
  function groupText(i) {
    const chem = C(), nb = chem.neighbors(st.mol, i).filter((j) => st.mol.atoms[j].el !== 'H');
    try { return chem.groupLabel(st.mol, i, nb.length ? nb[0] : null); } catch (e) { return st.mol.atoms[i].el; }
  }
  function chainTap(i, info) {
    const s = ch(), chem = C();
    if (s.showing || i == null || !st.mol || !st.mol.atoms[i]) return;
    info = info || {};
    if (s.phase === 'trace') {
      if (st.mol.atoms[i].el !== 'C') { chainSay('The main chain is carbons only. That ' + groupText(i) + ' is a group on the chain.', 'err'); flash([i], 'err'); renderChain(); return; }
      const last = s.path[s.path.length - 1];
      if (i === last) { if (!info.drag) chainUndo(); return; }
      if (s.path.indexOf(i) >= 0) { if (!info.drag) { chainSay('Already in your chain. Use Undo to back up.', 'err'); flash([i], 'err'); renderChain(); } return; }
      if (s.path.length && chem.bondBetween(st.mol, last, i) < 0) {
        if (info.drag) return;
        chainSay('Go bond by bond: tap a carbon bonded to the last one you picked.', 'err'); flash([i], 'err'); renderChain(); return;
      }
      s.path.push(i); s.flash = null;
      if ($('chain-msg').classList.contains('err')) chainSay('');
      renderChain();
      return;
    }
    if (s.phase === 'ring') {
      const res = s.res, ring = res.chain;
      if (s.ringC1 == null) {
        const err = ringFirstCheck(i);
        if (err) { chainSay(err, 'err'); flash([i], 'err'); renderChain(); return; }
        s.ringC1 = i; s.flash = null; chainSay('');
        renderChain();
        return;
      }
      if (i === s.ringC1) { s.ringC1 = null; renderChain(); return; }
      const order = chem.ringStart(st.mol, s.ringC1, i);
      if (!order.length || ring.indexOf(i) < 0) { chainSay('C2 has to be next to C1 in the ring.', 'err'); flash([i], 'err'); renderChain(); return; }
      const r = chem.checkNumbering(st.mol, order);
      if (r.ok) { acceptOrder(order, r.msg); return; }
      s.fails++;
      showWrongNumbers(order, [s.ringC1, i]);
      s.ringC1 = null;
      chainSay(r.msg, 'err');
      renderChain();
      return;
    }
    if (s.phase === 'direction') {
      if (i === s.path[0] || i === s.path[s.path.length - 1]) chainEnd(i);
    }
  }
  function ringFirstCheck(i) {
    const s = ch(), res = s.res, ring = res.chain;
    if (ring.indexOf(i) < 0) return 'Pick a carbon in the ring.';
    if (res.pcg) {
      if (res.pcg.atoms.indexOf(i) < 0) return 'C1 is the ring carbon with the ' + res.pcg.label + '.';
      return null;
    }
    if (res.unsat.length) {
      if (!res.unsat.some((u) => u.a === i || u.b === i)) return 'With a C=C in the ring, C1 and C2 are the two double-bond carbons.';
      return null;
    }
    if (res.substituents.length && !res.substituents.some((x) => x.at === i)) return 'Start at a carbon that has a branch.';
    return null;
  }
  function showWrongNumbers(order, atoms) {
    const s = ch(), dim = {};
    order.forEach((a, k) => { dim[a] = k + 1; });
    s.numbersDim = dim;
    flash(atoms, 'err');
    const tok = ++s.anim;
    setTimeout(() => { if (st.chain === s && s.anim === tok) { s.numbersDim = null; drawPad(); } }, 1500);
  }
  function chainCheck() {
    const s = ch(), chem = C();
    if (s.phase !== 'trace' || s.showing) return null;
    const r = chem.checkChain(st.mol, s.path);
    if (r.ok) {
      s.phase = 'direction'; s.flash = null;
      chainSay(r.msg, 'ok');
      renderChain();
      return r;
    }
    s.fails++;
    let msg = r.msg;
    if (r.code === 'shorter') {
      s.shorterFails++;
      if (s.shorterFails >= 2 && r.best && r.best.length) { flash([r.best[0]], 'hint'); msg += ' One end of a longer chain is marked.'; }
      else s.flash = null;
    } else if (r.code === 'end' || r.code === 'missing-pcg' || r.code === 'missing-unsat') flash(r.atoms, 'hint');
    else if (r.atoms && r.atoms.length) flash(r.atoms, 'err');
    else s.flash = null;
    chainSay(msg, 'err');
    renderChain();
    return r;
  }
  function chainEnd(e) {
    const s = ch(), chem = C();
    if (s.phase !== 'direction' || s.showing) return null;
    const p = s.path;
    if (e !== p[0] && e !== p[p.length - 1]) return null;
    const order = e === p[0] ? p.slice() : p.slice().reverse();
    const r = chem.checkNumbering(st.mol, order);
    if (r.ok) { acceptOrder(order, r.msg); return r; }
    s.fails++;
    showWrongNumbers(order, [e]);
    chainSay(r.msg, 'err');
    renderChain();
    return r;
  }
  function acceptOrder(order, msg) {
    const s = ch(), chem = C();
    const best = s.res.chain;
    if (sameArr(order, best)) { chem.setChainOrder(st.mol, null); if (st.src) chem.setChainOrder(st.src, null); }
    else if (chem.setChainOrder(st.mol, order)) { if (st.src) chem.setChainOrder(st.src, order); }
    s.order = st.mol.chainOrder ? order.slice() : null;
    s.phase = 'done'; s.tags = true; s.flash = null; s.numbersDim = null; s.path = []; s.ringC1 = null;
    relabel();
    chainSay(msg, 'ok');
    renderChain();
  }
  function chainUndo() {
    const s = ch();
    if (s.showing) return;
    if (s.phase === 'trace') { s.path.pop(); s.flash = null; }
    else if (s.phase === 'ring') s.ringC1 = null;
    else if (s.phase === 'direction') { s.phase = 'trace'; s.numbersDim = null; }
    chainSay('');
    renderChain();
  }
  function chainReset() {
    const s = ch();
    if (s.showing) return;
    s.path = []; s.ringC1 = null; s.flash = null; s.numbersDim = null;
    if (s.phase === 'direction') s.phase = 'trace';
    chainSay('');
    renderChain();
  }
  function chainAgain() {
    const s = ch(), chem = C();
    chem.setChainOrder(st.mol, null); if (st.src) chem.setChainOrder(st.src, null);
    s.order = null;
    relabel();
    chainStart();
  }
  function chainShowMe() {
    const s = ch();
    if (!s.res || s.showing) return;
    if (s.res.kind === 'single') { chainStart(); return; }
    if (s.phase !== 'trace' && s.phase !== 'ring' && s.phase !== 'direction') { s.phase = s.res.kind === 'chain' ? 'trace' : 'ring'; }
    const best = s.res.chain.slice(), tok = ++s.anim;
    s.showing = true; s.path = []; s.ringC1 = null; s.flash = null; s.numbersDim = null; s.showNumbers = {};
    chainSay('');
    const finish = () => {
      if (st.chain !== s || s.anim !== tok) return;
      s.showing = false; s.showNumbers = null;
      s.order = st.mol.chainOrder ? st.mol.chainOrder.slice() : null;
      s.phase = 'done'; s.tags = true; s.path = [];
      chainSay('Here\'s the main chain and the numbering.', 'ok');
      renderChain();
    };
    if (reduced()) { s.path = best.slice(); finish(); return; }
    let k = 0;
    const stepTrace = () => {
      if (st.chain !== s || s.anim !== tok) return;
      if (k < best.length) { s.path.push(best[k++]); renderChain(); setTimeout(stepTrace, 90); return; }
      k = 0; stepNum();
    };
    const stepNum = () => {
      if (st.chain !== s || s.anim !== tok) return;
      if (k < best.length) { s.showNumbers[best[k]] = k + 1; k++; drawPad(); setTimeout(stepNum, 60); return; }
      setTimeout(finish, 250);
    };
    renderChain();
    stepTrace();
  }
  function setTags(on) {
    ch().tags = !!on;
    renderChain();
  }

  // After the student's own (equivalent) numbering is accepted: everything that shows C numbers
  function relabel() {
    if (!st.mol) return;
    renderBondList();
    if (st.sel) {
      const chem = C();
      document.querySelectorAll('#bond-list .chip').forEach((b) => b.setAttribute('aria-pressed', String(+b.dataset.bond === st.sel.bond)));
      if (st.sel.info.rotatable) say('bond-msg', chem.bondLabel(st.mol, st.sel.front, st.sel.back) + ' selected. Look down it, then turn the back carbon.', 'ok');
      const keepMsg = $('prompt-msg').textContent;
      renderPrompt();
      if (keepMsg) say('prompt-msg', '');
      refresh(true);
    }
    buildTags();
  }

  /* ----- numbers on the 3D model ----- */
  let tagEls = [];
  function buildTags() {
    const layer = $('atom-tags');
    if (!layer) return;
    const s = ch(), chem = C();
    layer.innerHTML = ''; tagEls = [];
    const on = !!(s.tags && st.mol && s.res && (s.phase === 'done' || s.phase === 'skipped' || s.phase === 'intro'));
    layer.hidden = !on;
    if (!on) return;
    const order = mainOrder();
    st.mol.atoms.forEach((a, i) => {
      if (a.el !== 'C') return;
      const el = document.createElement('span');
      el.className = 'atom-tag' + (order.indexOf(i) < 0 ? ' branch' : '');
      el.textContent = chem.atomLabel(st.mol, i);
      layer.appendChild(el);
      tagEls.push({ i, el, r: (chem.EL.C.r) });
    });
  }
  function placeTags() {
    const layer = $('atom-tags');
    if (!layer || layer.hidden || !tagEls.length || !view || !st.local.length) return;
    const hide = st.overlayOn;
    let bs = null;
    try { bs = view.basis(); } catch (e) { bs = null; }
    tagEls.forEach((t) => {
      const p0 = st.local[t.i];
      if (!p0 || hide) { t.el.hidden = true; return; }
      const p = view.project(p0);
      let rad = 10;
      if (bs) { const q = view.project(V3.add(p0, V3.scale(bs.right, t.r))); rad = Math.hypot(q.x - p.x, q.y - p.y); }
      const behind = view.kind === 'gl' ? (p.z > 1 || p.z < -1) : false;
      if (behind || !isFinite(p.x) || !isFinite(p.y)) { t.el.hidden = true; return; }
      t.el.hidden = false;
      t.el.style.transform = 'translate(' + Math.round(p.x) + 'px,' + Math.round(p.y - rad - 10) + 'px) translate(-50%,-100%)';
    });
  }

  /* ------------------------------------------------------------------ */
  /* panel wiring                                                         */
  /* ------------------------------------------------------------------ */
  let sliderActive = false, lastHash = null;
  const sliderLast = {};
  function setInputMode(m) {
    ['type', 'draw', 'presets'].forEach((k) => {
      $('mode-' + k).setAttribute('aria-pressed', String(k === m));
      $('in-' + k).hidden = k !== m;
    });
    if (m === 'draw') ensureSketch();
  }
  function ensureSketch() {
    if (sketch || !window.NNSketch) { if (sketch) sketch.resize(); return; }
    sketch = NNSketch.attach($('sketch'), { mode: 'chain', bondLength: 44, picker: $('sk-picker') });
    sketch.onChange(onSketchChange);
    onSketchChange(sketch.getMolecule());
  }
  function onSketchChange(r) {
    $('sk-undo').disabled = !sketch || !sketch.canUndo();
    $('sk-clear').disabled = !sketch || sketch.isEmpty();
    if (!sketch || sketch.isEmpty()) { say('sk-msg', ''); $('sk-build').disabled = true; return; }
    if (r.ok) { sayHTML('sk-msg', subHTML(C().formula(r.mol)) + (r.mol.name ? ' (' + esc(r.mol.name) + ')' : '') + '. Ready to build.', 'ok'); $('sk-build').disabled = false; }
    else { say('sk-msg', r.error || r.message || 'Check the drawing.', 'err'); $('sk-build').disabled = true; }
  }
  function setSketchMode(m) {
    if (!sketch) return;
    sketch.setMode(m);
    $('sk-free').setAttribute('aria-pressed', String(m === 'free'));
    $('sk-chain').setAttribute('aria-pressed', String(m === 'chain'));
    $('sk-finish').hidden = m !== 'chain';
  }
  function fillPresets() {
    const sel = $('preset'), chem = C();
    if (!chem) return;
    const groups = [];
    chem.presets.forEach((p) => { if (!groups.includes(p.group)) groups.push(p.group); });
    sel.innerHTML = '<option value="" disabled selected>Choose a molecule…</option>' + groups.map((g) =>
      '<optgroup label="' + g.replace(/&/g, '&amp;') + '">' + chem.presets.filter((p) => p.group === g).map((p) =>
        '<option value="' + p.id + '">' + p.name + '</option>').join('') + '</optgroup>').join('');
  }
  function fromHash() {
    const h = decodeURIComponent((location.hash || '').slice(1));
    if (h === lastHash || '#' + h === lastHash) return false;
    if (!h) return false;
    const text = h.startsWith('smiles=') ? h.slice(7) : h;
    $('mol-input').value = text;
    return loadInput(text, { hash: false });
  }

  function wire() {
    $('mode-type').addEventListener('click', () => setInputMode('type'));
    $('mode-draw').addEventListener('click', () => setInputMode('draw'));
    $('mode-presets').addEventListener('click', () => setInputMode('presets'));
    $('mol-go').addEventListener('click', () => { const v = $('mol-input').value; loadInput(v); });
    $('mol-input').addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); loadInput($('mol-input').value); } });
    document.querySelectorAll('.hint .ex').forEach((b) => b.addEventListener('click', () => { $('mol-input').value = b.textContent; loadInput(b.textContent); }));
    $('preset').addEventListener('change', (e) => { if (e.target.value) loadInput(e.target.value, { msgId: 'mol-msg' }); });

    $('sk-free').addEventListener('click', () => setSketchMode('free'));
    $('sk-chain').addEventListener('click', () => setSketchMode('chain'));
    $('sk-finish').addEventListener('click', () => { if (sketch) sketch.setMode('chain'); });
    $('sk-undo').addEventListener('click', () => sketch && sketch.undo());
    $('sk-clear').addEventListener('click', () => sketch && sketch.clear());
    $('sk-build').addEventListener('click', () => {
      if (!sketch) return;
      const r = sketch.getMolecule();
      if (!r.ok) { say('sk-msg', r.error || 'Check the drawing.', 'err'); return; }
      if (loadMol(r.mol, { msgId: 'sk-msg' }) && window.matchMedia('(max-width: 860px)').matches) {
        document.querySelector('.viewer').scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' });
      }
    });

    $('look').addEventListener('click', () => lookDown(false));
    $('swap').addEventListener('click', () => {
      if (!st.sel) return;
      const wasOn = st.overlayOn;
      select(st.sel.bond, st.sel.back);
      if (wasOn) lookDown(false);
    });
    $('reset-view').addEventListener('click', () => { if (view && st.mol) view.reset(st.q0, reduced()); });

    // The sliders turn the back carbon by how far they move (relative), so the thumb never
    // jumps when the displayed pair of tied groups (two CH3 on one carbon) changes.
    ['dihedral', 'dihedral-m'].forEach((id) => {
      const sl = $(id);
      sl.addEventListener('input', () => {
        if (!st.sel) return;
        const v = +sl.value, last = sliderLast[id] != null ? sliderLast[id] : v;
        sliderActive = id; sliderLast[id] = v;
        const d = signed(v - last);
        if (d) rotateBy(st.sel.front, st.sel.back, d);
      });
      sl.addEventListener('change', () => { sliderActive = false; st.dirty = true; });
      sl.addEventListener('pointerdown', () => { st.dragging = true; sliderLast[id] = +sl.value; });
      sl.addEventListener('keydown', (e) => {
        if ((e.key === 'PageUp' || e.key === 'PageDown') && st.sel) {
          e.preventDefault();
          rotateBy(st.sel.front, st.sel.back, e.key === 'PageUp' ? 60 : -60);
        }
      });
    });
    window.addEventListener('pointerup', () => { if (st.dragging && sliderActive) { st.dragging = false; sliderActive = false; st.dirty = true; if (st.resultPending) writeResult(); } });
    document.querySelectorAll('.snap').forEach((b) => b.addEventListener('click', () => setDihedral(+b.dataset.deg)));
    $('flip').addEventListener('click', flipChair);

    $('prompt-check').addEventListener('click', checkPrompt);
    $('prompt-in').addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); checkPrompt(); } });
    $('prompt-next').addEventListener('click', () => { st.promptIdx++; renderPrompt(); });

    // main chain step
    $('chain-start').addEventListener('click', chainStart);
    $('chain-skip').addEventListener('click', chainSkip);
    $('chain-check').addEventListener('click', chainCheck);
    $('chain-undo').addEventListener('click', chainUndo);
    $('chain-reset').addEventListener('click', chainReset);
    $('chain-show').addEventListener('click', chainShowMe);
    $('chain-again').addEventListener('click', chainAgain);
    $('chain-open').addEventListener('click', () => { rememberSkip(false); chainStart(); });
    $('chain-tags').addEventListener('click', () => setTags(!ch().tags));
    $('chain-tags2').addEventListener('click', () => setTags(!ch().tags));

    // drag around the side-panel Newman circle
    if (N2()) {
      N2().attachDrag($('newman2d'), {
        getGeometry: () => (st.sel && st.sel.info.rotatable ? panelGeometry() : null),
        onStart: () => { st.dragging = true; $('newman2d').classList.add('dragging'); },
        onDelta: (d) => rotateBy(st.sel.front, st.sel.back, d),
        onEnd: () => { st.dragging = false; $('newman2d').classList.remove('dragging'); if (st.resultPending) writeResult(); }
      });
    }
    // redraw 2D canvases when their size changes
    const ro = new ResizeObserver(() => { st.dirty = true; });
    ro.observe($('newman2d')); ro.observe($('energy-plot'));
    window.addEventListener('hashchange', fromHash);
  }

  /* ------------------------------------------------------------------ */
  /* start                                                                */
  /* ------------------------------------------------------------------ */
  function loop() {
    requestAnimationFrame(loop);
    if (view) {
      checkOverlay();
      view.frame();
      placeTags();
    }
    if (st.dirty) refresh(false);
  }

  function start() {
    if (!window.NNChem || !window.NNGeom || !window.NNNewman2D) {
      $('mol-msg').textContent = 'Part of this tool did not load. Reload the page to try again.';
      $('mol-msg').className = 'msg err';
      return;
    }
    const host = $('canvas');
    try { view = hasWebGL() ? makeGL(host) : null; } catch (e) { console.warn('WebGL unavailable, using the flat renderer', e); view = null; host.innerHTML = ''; }
    if (!view) {
      st.flat = true;
      view = makeFlat(host);
      if (!view) {
        const d = document.createElement('div'); d.className = 'nogl';
        d.textContent = "Your browser can't draw 3D here. The Newman panel still works.";
        host.appendChild(d);
      }
    }
    if (view) setupViewerPointers();
    fillPresets();
    wire();
    if (!fromHash()) { $('mol-input').value = 'butane'; loadInput('butane', { hash: false }); }
    loop();
  }

  window.NNApp = {
    load: (input) => loadInput(String(input)),
    select: (bondIdx) => select(bondIdx),
    setDihedral: (deg) => { const ok = setDihedral(deg); refresh(true); return ok; },
    lookDown: (instant) => lookDown(instant),
    get state() { return st; },
    // main chain step (CHAIN_SPEC.md 6), for tests
    chain: {
      start: () => chainStart(),
      skip: () => chainSkip(),
      tap: (i) => chainTap(i, {}),
      check: () => chainCheck(),
      end: (i) => chainEnd(i),
      showMe: () => chainShowMe(),
      undo: () => chainUndo(),
      again: () => chainAgain(),
      state: () => st.chain,
      screenPos: (i) => (st.chain && st.chain.pad ? st.chain.pad.screenPos(i) : null)
    }
  };
  window.addEventListener('load', start);
})();

/*
 * rs/app.js: page controller for the R/S Assigner (rs.html).
 * Contract: rs/SPEC.md section 7. Only global: window.RSApp (test/debug handle).
 *
 * Uses NNChem / NNGeom / NNSketch (newman/, read-only), RSCip (engine), RSPresets, RSPaper, RSView3D.
 *
 * Deviations from SPEC.md (documented per the Newman convention):
 *  - Presets is the input tab that opens first (the default molecule is a preset); Type and Draw are one tap away.
 *  - The canvases sit in a .stage wrapper inside .viewer so the phone layout can give them one box.
 *  - Carbon chips and center names use NNChem.atomLabel, which now numbers carbons along the main chain
 *    (IUPAC-style), so chips, messages and the summary use the same numbers. Preset locants win for centers.
 *  - Step 2's "Yes" message lists the four groups in neighbor order, not priority order, so it doesn't give
 *    away step 3.
 *  - After a correct rank tap the tie-break table explains the pair that was just settled (k-1 vs k), not
 *    the next one, so it never reveals the next answer. Wrong taps still show (expected vs tapped).
 *  - Group labels that NNChem shortens (Bu, CR) are written out (CH2CH2CH2CH3, C(=CH2)CH3) when they fit.
 *  - A "Show me" button in step 5 traces 1 -> 2 -> 3 for the student (puts 4 in back first if needed).
 *  - Non-preset molecules get locants from NNChem.atomLabel and the NNChem name as the base, when it has one
 *    (else NNChem.mainChain's name). Carbons off the main chain are called "C8 (a branch carbon)" in messages and get
 *    no locant in the summary. Typed meso compounds use the IUPAC choice (R at the lower locant).
 *  - Step 5 asks "Which way does it turn?" before R/S; the 3D/paper direction word appears only after that answer.
 *    A wrong R/S answer gets a hint and a second try; the letter is shown only after a right answer or Show me.
 *  - Step 5 has chips for groups 1-3 (keyboard / phone tracing). On phones its controls move under the model (#dock).
 *  - Locked cards are inert. Group labels never end in a bare R ("CH(CH3)…"), and C(=O)OH reads COOH.
 *  - Centers the engine flags 'pseudo' or 'stereoDouble' are listed as special cases when the student taps None.
 */
(function () {
  'use strict';

  const $ = id => document.getElementById(id);
  const C = () => window.NNChem, G = () => window.NNGeom, R = () => window.RSCip, PR = () => window.RSPresets, PA = () => window.RSPaper;
  const reduced = () => !!(window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  const coarse = () => !!(window.matchMedia && window.matchMedia('(pointer: coarse)').matches);
  const narrow = () => !!(window.matchMedia && window.matchMedia('(max-width: 860px)').matches);
  const esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const subHTML = t => esc(t).replace(/([A-Za-z)\]])(\d+)/g, '$1<sub>$2</sub>');
  const cap = s => s.charAt(0).toUpperCase() + s.slice(1);
  const joinList = xs => (xs.length <= 1 ? xs.join('') : xs.slice(0, -1).join(', ') + ' and ' + xs[xs.length - 1]);
  const turnWord = t => (t === 'cw' ? 'clockwise' : 'counterclockwise');
  const flip = l => (l === 'R' ? 'S' : 'R');

  /* ------------------------------------------------------------------ */
  /* state                                                                */
  /* ------------------------------------------------------------------ */
  const st = {
    mol: null, coords: null, preset: null, want: null, mirrored: false, labels: {}, centers: [], found: new Set(),
    idx: -1, ranked: [], misses: {}, four: 'side', view: '3d', pose: 'dash', variant: 0, swaps: 0, newCount: 0,
    trace: [], turn: null, answers: {}, stage: 'none', names: {}, dist: null, slots: null,
    flash: {}, tb: null, arrowT: 0, traceSnap: null, findMisses: 0, warned: false, input: '', orientDirty: true, ready: false
  };
  let view = null, sketch = null, lastHash = null, flashSeq = 0;

  /* ------------------------------------------------------------------ */
  /* small UI helpers                                                     */
  /* ------------------------------------------------------------------ */
  function say(id, text, kind) { const e = $(id); if (!e) return; e.textContent = text || ''; e.className = 'msg' + (kind ? ' ' + kind : ''); }
  function sayHTML(id, html, kind) { const e = $(id); if (!e) return; e.innerHTML = html || ''; e.className = 'msg' + (kind ? ' ' + kind : ''); }
  function setCard(id, open) {
    const e = $(id);
    e.classList.toggle('locked', !open);
    e.inert = !open;              // locked cards stay out of the tab order
    if (id === 'c-trace') syncDock();
  }
  // phones: while step 5 is open, its controls sit under the model so the student can see both
  let dockTimer = null;
  function syncDock() {
    const body = $('trace-body'), dock = $('dock'), home = $('trace-home');
    if (!body || !dock || !home) return;
    const live = st.stage === 'orient' || st.stage === 'answer' || st.stage === 'done';
    const want = narrow() && live && !$('c-trace').classList.contains('locked');
    const docked = body.parentNode === dock;
    if (want && !docked) {
      dock.appendChild(body); dock.hidden = false; $('trace-moved').hidden = false;
      // bring the model and the controls into view, unless a tie-break explanation is open to read
      clearTimeout(dockTimer);
      if (st.ready && $('tiebreak').hidden) dockTimer = setTimeout(() => { if (body.parentNode === dock) showViewerOnPhone(); }, reduced() ? 0 : 500);
    } else if (!want && docked) {
      home.appendChild(body); dock.hidden = true; $('trace-moved').hidden = true;
    }
  }
  // desktop: bring a newly opened card into view inside the scrolling panel (phones keep their place)
  function reveal(id) {
    if (narrow() || !st.ready) return;
    const el = $(id); if (!el) return;
    try { el.scrollIntoView({ block: 'nearest', behavior: reduced() ? 'auto' : 'smooth' }); } catch (e) { /* old browsers */ }
  }
  function setDone(id, on) { $(id).classList.toggle('done', !!on); }
  function inputMode() { return ['type', 'draw', 'presets'].find(k => $('mode-' + k).getAttribute('aria-pressed') === 'true') || 'type'; }
  function msgId() { const m = inputMode(); return m === 'draw' ? 'sk-msg' : m === 'presets' ? 'preset-msg' : 'mol-msg'; }
  function setPending(on) { document.querySelector('.viewer').classList.toggle('pending', !!on); }

  /* ------------------------------------------------------------------ */
  /* names                                                                */
  /* ------------------------------------------------------------------ */
  const heavyNbrs = (i, from) => C().neighbors(st.mol, i).filter(j => st.mol.atoms[j].el !== 'H' && j !== from);
  const hCount = i => C().neighbors(st.mol, i).filter(j => st.mol.atoms[j].el === 'H').length;
  const bondOrder = (a, b) => { const k = C().bondBetween(st.mol, a, b); return k >= 0 ? (st.mol.bonds[k].order || 1) : 1; };
  const hPart = n => (n === 0 ? '' : n === 1 ? 'H' : 'H' + n);

  function condensed(s, from, seen) {
    seen.add(s);
    const a = st.mol.atoms[s];
    const base = a.el + hPart(hCount(s));
    const kids = heavyNbrs(s, from).filter(k => !seen.has(k)).map(k => {
      const o = bondOrder(s, k);
      return { t: (o === 2 ? '=' : o === 3 ? '≡' : '') + condensed(k, s, seen), o };
    });
    if (!kids.length) return base;
    if (kids.length > 1 && kids.every(k => k.t === kids[0].t)) return base + '(' + kids[0].t + ')' + kids.length;
    kids.sort((x, y) => (y.o - x.o) || (x.t.length - y.t.length));
    const last = kids.pop();
    return base + kids.map(k => '(' + k.t + ')').join('') + last.t;
  }
  // C(=O)OH -> COOH, C(=O)H -> CHO, C(=O)NH2 -> CONH2, so the same group reads the same everywhere
  const tidy = t => t.replace(/C\(=O\)OH(?![a-z\d])/g, 'COOH').replace(/C\(=O\)NH2/g, 'CONH2').replace(/C\(=O\)H(?![a-z\d])/g, 'CHO');
  function niceLabel(s, from) {
    if (st.mol.atoms[s].el === 'H') return 'H';
    let lab;
    try { lab = C().groupLabel(st.mol, s, from); } catch (e) { lab = st.mol.atoms[s].el; }
    if (/^(Et|Pr|iPr|Bu|iBu|sBu|tBu)$/.test(lab) || /[A-Z][a-z]?H?\d*R\d*$/.test(lab)) {
      try {
        const t = tidy(condensed(s, from, new Set([from])));
        if (t.length <= 20) lab = t;
        else {
          // first atom and its end groups, then "…" for the rest of the chain: CH(OH)…
          // (never a bare R here: in an R/S tool it reads as the descriptor)
          const a = st.mol.atoms[s];
          const kids = heavyNbrs(s, from);
          const ends = kids.filter(k => heavyNbrs(k, s).length === 0).map(k => (bondOrder(s, k) === 2 ? '=' : '') + st.mol.atoms[k].el + hPart(hCount(k)));
          const rest = kids.length - ends.length;
          lab = a.el + hPart(hCount(s)) + ends.map(e => '(' + e + ')').join('') + (rest ? '…' : '');
        }
      } catch (e) { /* keep NNChem's label */ }
    }
    return tidy(lab).replace(/R\d*$/, '…');
  }
  function centerName(i) {
    const p = st.preset;
    if (p && p.locants && p.locants[i] != null) return 'C' + p.locants[i];
    try { return C().atomLabel(st.mol, i); } catch (e) { return 'C' + (i + 1); }
  }
  function centerNum(i) {
    const p = st.preset;
    if (p && p.locants && p.locants[i] != null) return p.locants[i];
    const m = /^C(\d+)$/.exec(centerName(i));
    return m ? +m[1] : 1000 + i;
  }
  function carbonName(i) { return centerName(i); }
  // for messages: a carbon off the main chain gets a number after the chain (C8 on a heptane), so say it is a branch
  function onMainChain(i) {
    if (st.preset && st.preset.locants && st.preset.locants[i] != null) return true;
    if (st.chainSet === undefined) {
      st.chainSet = null;
      try { const mc = C().mainChain(st.mol); if (mc && mc.chain && mc.chain.length) st.chainSet = new Set(mc.chain); } catch (e) { /* no chain */ }
    }
    return !st.chainSet || st.chainSet.has(i);
  }
  function cText(i) { return centerName(i) + (onMainChain(i) ? '' : ' (a branch carbon)'); }

  function groupNames(c) {
    const out = {};
    const ligs = c.ranking.ligands;
    ligs.forEach(l => { out[l.atom] = (l.ring ? 'ring ' : '') + niceLabel(l.atom, c.atom); });
    for (let i = 0; i < ligs.length; i++) {
      for (let j = i + 1; j < ligs.length; j++) {
        const a = ligs[i].atom, b = ligs[j].atom;
        if (out[a] !== out[b]) continue;
        const cmp = R().compareLigands(st.mol, c.atom, a, b);
        if (!cmp.sign) continue;
        const hi = cmp.sign > 0 ? a : b, lo = hi === a ? b : a, side = hi === a ? 'a' : 'b';
        out[hi] += ' (toward ' + feature(cmp, side) + ')';
        out[lo] += ' (other side)';
      }
    }
    return out;
  }
  function feature(cmp, side) {
    const box = cmp.boxes && cmp.boxes[side];
    const dup = box && box.atoms.find(e => e.dup);
    if (dup) return dup.el === 'O' ? 'C=O' : dup.el === 'C' ? 'C=C' : 'C=' + dup.el;
    const d = cmp.deciding && cmp.deciding[side];
    if (d && d.el !== 'C' && d.el !== 'H' && d.el !== '0') return 'the ' + d.el;
    return 'the branch';
  }
  const gname = a => st.names[a] || '?';
  const gHTML = a => subHTML(gname(a));

  // msg with {C}, {X}, {Y} placeholders -> HTML
  function fill(msg, map) {
    let h = esc(msg);
    Object.keys(map).forEach(k => { h = h.split('{' + k + '}').join(map[k]); });
    return h;
  }

  /* ------------------------------------------------------------------ */
  /* loading                                                              */
  /* ------------------------------------------------------------------ */
  function loadText(text, opts) {
    opts = Object.assign({ msgId: 'mol-msg' }, opts || {});
    text = String(text || '').trim();
    $('mol-choices').hidden = true; setPending(false);
    if (!text) { say(opts.msgId, 'Type a name, a condensed formula, or SMILES first.', 'err'); return false; }
    if (!C() || !R()) { say(opts.msgId, 'Part of this tool did not load. Reload the page.', 'err'); return false; }
    if (/\[\d+[A-Z]|\[D\]/.test(text)) { say(opts.msgId, "This tool doesn't handle isotopes like D (²H). Try the molecule without the isotope label.", 'err'); return false; }
    const lk = PR().lookup(text);
    if (lk) return loadPreset(lk.preset, Object.assign({}, opts, { want: lk.want, input: text }));
    const t2 = text.replace(/[\u2212\u2013\u2014]/g, '-');
    const m = /@/.test(t2) ? null : t2.match(PR().DESC);
    const want = m ? '(' + m[1].split(',').map(s => s.trim().toUpperCase()).join(',') + ')' : null;
    let r;
    try { r = C().parse(text); } catch (e) { console.warn('NNChem.parse failed', e); r = { ok: false }; }
    if (!r.ok && m) { try { r = C().parse(t2.slice(m[0].length).trim()); } catch (e) { r = { ok: false }; } }
    if (r.ok) return loadMol(r.mol, Object.assign({}, opts, { want, input: text }));
    if (r.choices) {
      say(opts.msgId, r.message || 'That fits more than one molecule. Which one?');
      const box = $('mol-choices');
      box.innerHTML = '';
      r.choices.forEach(ch => {
        const b = document.createElement('button');
        b.type = 'button'; b.className = 'chip'; b.textContent = ch.label;
        b.addEventListener('click', () => { $('mol-input').value = ch.value; loadText(ch.value); });
        box.appendChild(b);
      });
      box.hidden = false; setPending(true);
      return false;
    }
    say(opts.msgId, r.error || 'I could not read that. Try a name, a condensed formula, or SMILES.', 'err');
    return false;
  }

  function loadPreset(p, opts) {
    opts = opts || {};
    const r = C().parse(p.smiles);
    if (!r.ok) { say(opts.msgId || msgId(), 'That preset did not load.', 'err'); return false; }
    return loadMol(r.mol, Object.assign({}, opts, { preset: p }));
  }

  function descriptorOf(labels, preset) {
    const atoms = Object.keys(labels).map(Number);
    if (!atoms.length) return '';
    if (atoms.length === 1) return '(' + labels[atoms[0]] + ')';
    const loc = preset && preset.locants ? preset.locants : null;
    if (!loc || atoms.some(a => loc[a] == null)) return null;
    atoms.sort((x, y) => loc[x] - loc[y]);
    return '(' + atoms.map(a => loc[a] + labels[a]).join(',') + ')';
  }
  const invertDesc = d => d.replace(/[RS]/g, ch => (ch === 'R' ? 'S' : 'R'));
  const sameDesc = (a, b, n) => (n === 1 ? a.replace(/\d/g, '') === b.replace(/\d/g, '') : a === b);

  function build(mol) {
    const b = G().build3D(mol);
    const centers = R().findStereocenters(b.mol, b.coords);
    const labels = {}; centers.forEach(c => { labels[c.atom] = c.label; });
    return { mol: b.mol, coords: b.coords, centers, labels };
  }

  function loadMol(mol, opts) {
    opts = opts || {};
    const mid = opts.msgId || msgId();
    let b;
    try { b = build(mol); } catch (e) {
      console.warn('build failed', e);
      say(mid, 'I could not build a 3D model of that molecule. Try a smaller or simpler one.', 'err');
      return false;
    }
    let preset = opts.preset || null, mirrored = false;
    if (opts.want) {
      const n = b.centers.length;
      const fits = (bb, p) => {
        const d = descriptorOf(bb.labels, p);
        if (d == null) return null;
        if (sameDesc(d, opts.want, n)) return 'same';
        if (sameDesc(invertDesc(d), opts.want, n)) return 'mirror';
        return null;
      };
      let f = n ? fits(b, preset) : null;
      if (!f && preset && n > 1) {
        // try the other presets with the same base (meso vs chiral diastereomers)
        const others = PR().list.filter(p => p.base === preset.base && p.id !== preset.id);
        for (const p of others) {
          const bb = build(C().parse(p.smiles).mol);
          const ff = fits(bb, p);
          if (ff) { preset = p; b = bb; f = ff; break; }
        }
      }
      if (!f) {
        say(mid, n === 0 ? 'That molecule has no stereocenter, so it has no R or S.'
          : preset ? "That's a diastereomer of the preset. I can only build the preset or its mirror image from a name. Type SMILES with @ and @@ instead."
            : "I can't build that exact stereoisomer from a name. Type SMILES with @ and @@ instead.", 'err');
        if (n !== 0) return false;
      }
      if (f === 'mirror') {
        b.coords = R().mirror(b.coords);
        b.centers = R().findStereocenters(b.mol, b.coords);
        b.labels = {}; b.centers.forEach(c => { b.labels[c.atom] = c.label; });
        mirrored = true;
      }
    }
    st.mol = b.mol; st.coords = b.coords; st.preset = preset; st.mirrored = mirrored; st.want = opts.want || null; st.chainSet = undefined;
    st.input = opts.input || (preset ? preset.name : (mol.input || ''));
    st.src = mol;
    computeCenters();
    if (view) view.setMolecule(st.mol, st.coords);

    // message
    let what;
    if (preset) {
      const d = R().describe(st.mol, st.labels, { locants: preset.locants, base: preset.base });
      what = esc(mirrored && !d.meso && d.text ? d.text : preset.name);
    } else {
      const fml = subHTML(C().formula(st.mol));
      const nm = cleanName(st.mol.name || mol.name);
      what = nm ? esc(nm) + ' (' + fml + ')' : fml;
    }
    if (opts.msgId !== false) sayHTML(mid, 'Built ' + what + '.', 'ok');
    $('stereo-note').hidden = !st.centers.some(c => !c.given);
    $('mirror').disabled = false;
    if (preset) $('preset').value = preset.id; else $('preset').value = '';

    if (opts.hash !== false) {
      let h = '';
      if (preset && !mirrored) h = '#' + preset.id;
      else if (preset || (mol.source !== 'sketch' && st.input)) h = '#' + (preset ? encodeURIComponent(st.input) : 'smiles=' + encodeURIComponent(mol.input || st.input));
      try { history.replaceState(null, '', location.pathname + location.search + h); } catch (e) { /* file:// */ }
      lastHash = h;
    }
    renderLegend();
    resetAll();
    return true;
  }
  function cleanName(n) {
    if (!n) return '';
    return String(n).replace(/^meso-/i, '').replace(/^\((?:\d?[RS])(?:,\s*\d?[RS])*\)-/i, '');
  }

  function computeCenters() {
    const cs = R().findStereocenters(st.mol, st.coords);
    cs.sort((a, b) => (centerNum(a.atom) - centerNum(b.atom)) || (a.atom - b.atom));
    st.centers = cs;
    st.labels = {}; cs.forEach(c => { st.labels[c.atom] = c.label; });
    st.special = [];
    for (let i = 0; i < st.mol.nHeavy; i++) {
      if (st.mol.atoms[i].el !== 'C') continue;
      try { const k = R().classifyAtom(st.mol, i, st.coords).code; if (k === 'pseudo' || k === 'stereoDouble') st.special.push(i); } catch (e) { /* skip */ }
    }
  }

  /* ------------------------------------------------------------------ */
  /* resets                                                               */
  /* ------------------------------------------------------------------ */
  function resetAll() {
    st.found = new Set(); st.findMisses = 0;
    st.stage = 'find';
    resetCenterSteps();
    st.idx = -1; st.answers = {};
    setCard('c-find', true); setDone('c-find', false);
    ['c-rank', 'c-orient', 'c-trace', 'c-name'].forEach(id => { setCard(id, false); setDone(id, false); });
    say('find-msg', '');
    $('summary-text').textContent = 'Finish each stereocenter first.';
    $('summary-text').classList.add('wait');
    $('summary-note').textContent = '';
    $('result').hidden = true;
    $('view-paper').disabled = true;
    if (st.view !== '3d') setView('3d');
    renderAtomList();
    $('rank-q').textContent = 'Tap the groups from highest priority (1) to lowest (4).';
    $('group-list').innerHTML = '';
    $('orient-state').textContent = 'Rank the groups first.';
    if (view) view.setPickable(null);
    updateStyle();
  }
  function resetCenterSteps() {
    st.ranked = []; st.misses = {}; st.trace = []; st.turn = null; st.swaps = 0; st.pose = 'dash'; st.variant = 0; st.newCount = 0;
    st.tb = null; st.arrowT = 0; st.traceSnap = null; st.arrowHidden = false; st.paperWrong = null;
    say('rank-msg', ''); say('trace-msg', '');
    $('tiebreak').hidden = true; $('tiebreak').open = false;
    $('rank-show').className = 'btn-ghost';
    st.turnKnown = false; st.ansMissed = false;
    $('turn').hidden = true; $('turn-row').hidden = true; $('answer-row').hidden = true; $('trace-tools').hidden = false; $('next-center').hidden = true;
    $('trace-chips').hidden = false; $('trace-chips').innerHTML = '';
    resetPickRows();
    $('swap').hidden = true;
  }

  /* ------------------------------------------------------------------ */
  /* step 2: find                                                         */
  /* ------------------------------------------------------------------ */
  function renderAtomList() {
    const box = $('atom-list');
    box.innerHTML = '';
    if (!st.mol) return;
    const carbons = [];
    for (let i = 0; i < st.mol.nHeavy; i++) if (st.mol.atoms[i].el === 'C') carbons.push(i);
    carbons.sort((a, b) => (centerNum(a) - centerNum(b)) || (a - b));
    carbons.forEach(i => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip'; b.dataset.atom = i;
      b.textContent = carbonName(i);
      b.setAttribute('aria-pressed', String(st.found.has(i)));
      b.addEventListener('click', () => tapAtom(i));
      box.appendChild(b);
    });
  }
  function syncAtomChips() {
    document.querySelectorAll('#atom-list .chip').forEach(b => b.setAttribute('aria-pressed', String(st.found.has(+b.dataset.atom))));
  }
  function ligandAtoms(center, lig) {
    // the atoms of a group: everything reachable from lig without passing through center
    const out = [lig], seen = new Set([center, lig]), stack = [lig];
    while (stack.length) {
      const x = stack.pop();
      C().neighbors(st.mol, x).forEach(y => { if (!seen.has(y)) { seen.add(y); out.push(y); stack.push(y); } });
    }
    return out;
  }
  function flash(atoms, style, ms) {
    const key = ++flashSeq;
    atoms.forEach(a => { st.flash[a] = { style, key }; });
    updateStyle();
    setTimeout(() => { atoms.forEach(a => { if (st.flash[a] && st.flash[a].key === key) delete st.flash[a]; }); updateStyle(); }, ms);
  }

  function tapAtom(i) {
    if (!st.mol || i == null) return;
    if (st.stage !== 'find') {
      if (st.stage === 'rank' || st.stage === 'orient' || st.stage === 'answer') return tapGroup(groupFor(i));
      return;
    }
    if (st.found.has(i)) {
      st.found.delete(i); syncAtomChips(); updateStyle();
      say('find-msg', cText(i) + ' unmarked.');
      return;
    }
    const cls = R().classifyAtom(st.mol, i, st.coords);
    const name = esc(cText(i));
    if (cls.code === 'ok') {
      st.found.add(i);
      const nb = C().neighbors(st.mol, i);
      const list = nb.map(l => subHTML((R().rankSubstituents(st.mol, i).ligands.find(x => x.atom === l).ring ? 'ring ' : '') + niceLabel(l, i)));
      sayHTML('find-msg', 'Yes. ' + name + ' has four different groups: ' + joinList(list) + '.', 'ok');
      syncAtomChips(); updateStyle();
      return;
    }
    let msg = cls.msg;
    if (cls.code === 'sameGroups' && cls.pair) {
      const A = niceLabel(cls.pair[0], i);
      sayHTML('find-msg', 'This carbon has two identical groups (' + subHTML(A) + ' and ' + subHTML(A) + '). A stereocenter needs four different groups.', 'err');
    } else {
      sayHTML('find-msg', fill(msg, { C: name }), 'err');
    }
    flash([i], 'wrong', 900);
    if ((cls.code === 'sameGroups' || cls.code === 'ringSame') && cls.pair) {
      const A = ligandAtoms(i, cls.pair[0]), B = ligandAtoms(i, cls.pair[1]);
      st.findPair = { a: A, b: B, until: Date.now() + 1500 };
      updateStyle();
      setTimeout(() => { if (st.findPair && Date.now() >= st.findPair.until - 5) { st.findPair = null; updateStyle(); } }, 1500);
    }
  }

  function centerList(atoms) { return joinList(atoms.map(a => cText(a))); }
  function findDone() {
    if (st.stage !== 'find') return;
    const n = st.centers.length;
    const got = st.centers.filter(c => st.found.has(c.atom)).length;
    if (n === 0) {
      if (st.found.size === 0) return findNone();
      say('find-msg', 'None of these has four different groups. Unmark them.', 'err');
      return;
    }
    if (got === n) {
      say('find-msg', (n === 1 ? "That's the only one: " + cText(st.centers[0].atom) + '.' : "That's all of them: " + centerList(st.centers.map(c => c.atom)) + '.'), 'ok');
      setDone('c-find', true);
      $('view-paper').disabled = false;
      startCenter(0);
      return;
    }
    st.findMisses++;
    let m = 'You have ' + got + ' of ' + n + '. Look again for a carbon with four different groups.';
    if (st.findMisses >= 2 && st.centers.some(c => !st.found.has(c.atom) && (st.mol.rings || []).some(r => r.indexOf(c.atom) >= 0))) m += ' Check the ring carbons.';
    say('find-msg', m, 'err');
  }
  function findNone() {
    if (st.stage !== 'find') return;
    if (st.centers.length) { say('find-msg', 'There is at least one. Look for a carbon with four different groups.', 'err'); return; }
    const note = st.preset && st.preset.note ? ' ' + st.preset.note : '';
    const sp = st.special || [];
    if (sp.length) {
      say('find-msg', 'Right, none that this tool can label. ' + joinList(sp.map(cText)) + (sp.length > 1 ? ' are special cases' : ' is a special case') +
        ': tap ' + (sp.length > 1 ? 'one' : 'it') + ' to see why.', 'ok');
    } else say('find-msg', 'Right. No carbon here has four different groups, so there is no R or S.' + note, 'ok');
    setDone('c-find', true);
    st.stage = 'finished';
    setCard('c-name', true); setDone('c-name', true);
    $('summary-text').textContent = sp.length ? 'No stereocenters this tool can label.' : 'No stereocenters.';
    $('summary-text').classList.remove('wait');
    $('summary-note').textContent = sp.length ? 'Some molecules need extra CIP rules (E/Z or r/s) that this tool skips.'
      : st.preset && st.preset.note ? st.preset.note : 'A molecule needs at least one carbon with four different groups to get R or S.';
  }
  function findShow() {
    if (st.stage !== 'find') return;
    if (!st.centers.length) { findNone(); return; }
    st.centers.forEach(c => st.found.add(c.atom));
    syncAtomChips(); updateStyle();
    const n = st.centers.length;
    say('find-msg', (n === 1 ? 'Here it is: ' + cText(st.centers[0].atom) + '. It has four different groups.'
      : 'Here they are: ' + centerList(st.centers.map(c => c.atom)) + '. Each has four different groups.') + " Tap That's all of them to go on.", 'ok');
  }

  /* ------------------------------------------------------------------ */
  /* step 3: rank                                                         */
  /* ------------------------------------------------------------------ */
  const cur = () => st.centers[st.idx] || null;
  const order = () => cur().ranking.order;

  function startCenter(k) {
    st.idx = k;
    const c = cur();
    if (!c) return;
    resetCenterSteps();
    st.stage = 'rank';
    st.names = groupNames(c);
    // graph distances from each ligand, not passing through the center
    st.dist = {};
    c.ranking.ligands.forEach(l => {
      const d = {}; d[l.atom] = 0;
      const q = [l.atom], seen = new Set([c.atom, l.atom]);
      while (q.length) { const x = q.shift(); C().neighbors(st.mol, x).forEach(y => { if (!seen.has(y)) { seen.add(y); d[y] = d[x] + 1; q.push(y); } }); }
      st.dist[l.atom] = d;
    });
    computeSlots();
    renderGroups();
    const nm = cText(c.atom);
    const multi = st.centers.length > 1 ? ' (' + (k + 1) + ' of ' + st.centers.length + ')' : '';
    $('rank-q').textContent = 'Stereocenter ' + nm + multi + '. Tap the groups from highest priority (1) to lowest (4).';
    setCard('c-rank', true); setDone('c-rank', false);
    ['c-orient', 'c-trace'].forEach(id => { setCard(id, false); setDone(id, false); });
    setCard('c-name', false);
    $('orient-state').textContent = 'Rank the groups first.';
    $('result').hidden = false;
    $('r-label').textContent = '?';
    $('r-center').textContent = nm;
    $('r-prio').textContent = 'not ranked yet';
    $('r-four').textContent = '–';
    $('r-summary').textContent = st.centers.length > 1 ? Object.keys(st.answers).length + ' of ' + st.centers.length + ' stereocenters done' : '1 stereocenter';
    if (view) {
      const pick = []; for (let i = 0; i < st.mol.atoms.length; i++) if (i !== c.atom) pick.push(i);
      view.setPickable(pick);
    }
    updateStyle(); drawPaper(); updateOrientation(true);
    reveal('c-rank');
  }

  function renderGroups() {
    const c = cur(), box = $('group-list');
    box.innerHTML = '';
    c.ranking.ligands.forEach(l => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip'; b.dataset.atom = l.atom;
      b.innerHTML = '<span class="gl">' + gHTML(l.atom) + '</span><span class="badge"></span>';
      b.addEventListener('click', () => tapGroup(l.atom));
      box.appendChild(b);
    });
    syncGroupChips();
  }
  function syncGroupChips() {
    document.querySelectorAll('#group-list .chip').forEach(b => {
      const a = +b.dataset.atom, k = st.ranked.indexOf(a);
      b.querySelector('.badge').textContent = k >= 0 ? String(k + 1) : '';
      b.classList.toggle('ranked', k >= 0);
      b.setAttribute('aria-label', gname(a) + (k >= 0 ? ', priority ' + (k + 1) : ''));
    });
  }

  function groupFor(atom) {
    const c = cur();
    if (!c || atom == null || atom === c.atom) return null;
    const ligs = c.ranking.ligands.map(l => l.atom);
    if (ligs.indexOf(atom) >= 0) return atom;
    let best = [], bd = Infinity;
    ligs.forEach(l => { const d = st.dist[l][atom]; if (d == null) return; if (d < bd) { bd = d; best = [l]; } else if (d === bd) best.push(l); });
    if (best.length <= 1) return best[0] != null ? best[0] : null;
    if (!view) return best[0];
    const p = view.project(atom);
    best.sort((x, y) => { const a = view.project(x), b = view.project(y); return Math.hypot(a.x - p.x, a.y - p.y) - Math.hypot(b.x - p.x, b.y - p.y); });
    return best[0];
  }

  function tapGroup(lig) {
    if (lig == null || !cur()) return;
    if (st.stage === 'rank') return rankTap(lig);
    if (st.stage === 'orient') {
      if (st.view === '3d' && view && view.busy()) { say('trace-msg', 'Wait for the model to stop turning, then trace.'); return; }
      if (traceAllowed()) return traceTap(lig);
      say('trace-msg', st.view === 'paper' ? 'Get group 4 onto the dash (or the wedge) first.' : 'Turn the model so group 4 points away from you (or toward you) first.', 'err');
      return;
    }
    if (st.stage === 'answer') say('trace-msg', st.turnKnown ? 'Now pick R or S.' : 'First say which way the arrow turns.', 'err');
  }

  function rankTap(X) {
    const k = st.ranked.length + 1;
    if (k > 4) return;
    if (st.ranked.indexOf(X) >= 0) { sayHTML('rank-msg', gHTML(X) + ' already has number ' + (st.ranked.indexOf(X) + 1) + '. Next is ' + k + '.', 'err'); return; }
    const Y = order()[k - 1];
    st.tb = null;
    if (X === Y) { accept(k, false); return; }
    st.misses[k] = (st.misses[k] || 0) + 1;
    const d = R().diagnose(st.mol, cur().atom, X, Y, k);
    sayHTML('rank-msg', fill(d.msg, { C: esc(centerName(cur().atom)), X: gHTML(X), Y: gHTML(Y) }), 'err');
    showTable(d.comparison, Y, X, false);
    const chip = document.querySelector('#group-list .chip[data-atom="' + X + '"]');
    if (chip) {
      chip.classList.remove('shake'); void chip.offsetWidth; chip.classList.add('shake', 'wrong');
      setTimeout(() => chip.classList.remove('shake', 'wrong'), 900);
    }
    st.paperWrong = X; drawPaper();
    setTimeout(() => { if (st.paperWrong === X) { st.paperWrong = null; drawPaper(); } }, 900);
    flash([X], 'wrong', 900);
    if (st.misses[k] >= 2) $('rank-show').className = 'btn';
  }

  function accept(k, shown) {
    const ord = order();
    const Y = ord[k - 1];
    st.ranked.push(Y);
    let html = shown ? 'Number ' + k + ' is ' + gHTML(Y) + '.' : (k === 1 ? 'Yes, ' + gHTML(Y) + ' is 1.' : 'Yes, ' + k + '.');
    if (k === 3) {
      st.ranked.push(ord[3]);
      html += ' That leaves ' + gHTML(ord[3]) + ' as 4.';
    }
    $('rank-show').className = 'btn-ghost';
    sayHTML('rank-msg', html, 'ok');
    // explain the pair that was just settled (k-1 vs k), and on the last step also 3 vs 4
    const steps = cur().ranking.steps;
    let shownTb = false;
    const showStep = j => { const s = steps[j]; if (s && s.sphere >= 2) { showTable(s, ord[j], ord[j + 1], true); shownTb = true; } };
    if (k === 3) { showStep(2); if (!shownTb) showStep(1); } else if (k >= 2) showStep(k - 2);
    if (!shownTb) { $('tiebreak').hidden = true; st.tb = null; }
    syncGroupChips(); updateStyle(); drawPaper();
    if (st.ranked.length === 4) finishRank();
  }

  function finishRank() {
    const ord = order();
    $('r-prio').innerHTML = ord.map(a => gHTML(a)).join(' &gt; ');
    setDone('c-rank', true);
    st.stage = 'orient';
    renderTraceChips();
    setCard('c-orient', true);
    updateOrientation(true);
    drawPaper();
    reveal('c-trace');
    if (st.view === '3d' && narrow()) { /* keep the panel where it is on phones */ }
  }

  function rankUndo() {
    if (st.stage !== 'rank' && st.stage !== 'orient') return;
    if (!st.ranked.length) return;
    if (st.ranked.length === 4) { st.ranked.pop(); st.ranked.pop(); } else st.ranked.pop();
    if (st.stage === 'orient') {
      st.stage = 'rank'; clearTrace(); $('trace-chips').innerHTML = '';
      setCard('c-orient', false); setCard('c-trace', false); setDone('c-rank', false); setDone('c-orient', false);
      $('r-prio').textContent = 'not ranked yet';
      $('orient-state').textContent = 'Rank the groups first.';
    }
    st.tb = null; $('tiebreak').hidden = true;
    say('rank-msg', st.ranked.length ? 'Undone. Next is ' + (st.ranked.length + 1) + '.' : 'Undone. Start with 1.');
    syncGroupChips(); updateStyle(); drawPaper();
  }
  function rankShow() {
    if (st.stage !== 'rank') return;
    accept(st.ranked.length + 1, true);
    if (st.stage === 'rank') return;
  }

  /* ---- tie-break table ---- */
  function boxHTML(box, diffPos) {
    if (!box) return '';
    if (box.from == null) {
      const e = box.atoms[0];
      return '<span class="' + (diffPos === 0 ? 'd' : '') + '">' + esc(e.el) + '</span>';
    }
    return '(' + box.atoms.map((e, j) => {
      const cls = [j === diffPos ? 'd' : '', e.dup ? 'dup' : ''].filter(Boolean).join(' ');
      const title = e.dup ? ' title="duplicate (counts the multiple bond again)"' : e.el === '0' ? ' title="nothing there (phantom atom)"' : '';
      return '<span' + (cls ? ' class="' + cls + '"' : '') + title + '>' + esc(e.el) + '</span>';
    }).join(',') + ')';
  }
  function showTable(cmp, yAtom, xAtom, settled) {
    if (!cmp || !cmp.rows) return;
    const rows = cmp.rows.map(row => {
      const keep = (bx, k) => !bx.terminal || (row.diff && row.diff.box === k);
      const cell = side => row[side].map((bx, k) => (keep(bx, k) ? boxHTML(bx, row.diff && row.diff.box === k ? row.diff.pos : -1) : null)).filter(x => x != null).join(' ') || '<span class="tie">–</span>';
      return '<tr' + (row.diff ? ' class="diff"' : '') + '><td>' + row.sphere + '</td><td>' + cell('a') + '</td><td>' + cell('b') + '</td></tr>';
    }).join('');
    const ya = gname(yAtom), xa = gname(xAtom);
    const steps = R().explain ? R().explain(cmp, { a: ya, b: xa }) : [];
    $('tb-table').innerHTML = '<table><thead><tr><th>Sphere</th><th class="ga">' + subHTML(ya) + '</th><th class="gb">' + subHTML(xa) + '</th></tr></thead><tbody>' + rows + '</tbody></table>' +
      (steps.length ? '<ul class="tb-steps">' + steps.map(s => '<li>' + subHTML(s) + '</li>').join('') + '</ul>' : '') +
      '<p class="tb-key">Sphere 1 is the atom attached to ' + esc(centerName(cur().atom)) + '. Underlined italic = a copy (from a double or triple bond, or where a ring closes back on itself). 0 = nothing there.</p>';
    $('tiebreak').hidden = false;
    $('tiebreak').open = true;
    // one <span> inside the flex summary, so <sub> stays a subscript
    $('tiebreak').querySelector('summary').innerHTML = '<span>' + (settled ? 'Why ' + subHTML(ya) + ' outranks ' + subHTML(xa) : 'How the tie is broken') + '</span>';
    // 3D highlight of the deciding sets (teal = the winner, pink = the loser) until the next tap
    const atomsOf = box => (box ? [box.from].concat(box.atoms.filter(e => e.atom != null && !e.dup).map(e => e.atom)).filter(a => a != null) : []);
    if (cmp.boxes && cmp.sphere >= 2) st.tb = { a: atomsOf(cmp.boxes.a), b: atomsOf(cmp.boxes.b) };
    else st.tb = { a: [yAtom], b: [xAtom] };
    updateStyle();
  }

  /* ------------------------------------------------------------------ */
  /* step 4: orient                                                       */
  /* ------------------------------------------------------------------ */
  function computeSlots() {
    const c = cur(); if (!c) return;
    const ord = c.ranking.order;
    const vecs = ord.map(a => st.coords[a].map((x, k) => x - st.coords[c.atom][k]));
    st.slots = PA().layout(vecs, st.pose, st.variant);
    if (st.swaps) {
      const j = st.slots.indexOf('dash'); const s4 = st.slots[3]; st.slots[3] = 'dash'; st.slots[j] = s4;
    }
  }
  // reads on its own (the phone result card hides the 'Group 4' key)
  function fourText(four) {
    if (four === 'back') return '4 away from you';
    if (four === 'front') return '4 toward you';
    return '4 sideways';
  }
  function updateOrientation(force) {
    const c = cur();
    if (!c || (st.stage !== 'orient' && st.stage !== 'answer' && st.stage !== 'done')) return;
    const lab4 = gHTML(order()[3]);
    let four, html;
    if (st.view === 'paper') {
      four = st.pose === 'dash' ? 'back' : st.pose === 'wedge' ? 'front' : (st.swaps ? 'back' : 'side');
      if (st.pose === 'plane' && st.swaps) html = 'I swapped 4 with the dash group. That makes the mirror image, so remember to flip at the end.';
      else if (four === 'back') html = 'Group 4 (' + lab4 + ') is on the dash, pointing away from you. Go to step 5.';
      else if (four === 'front') html = 'Group 4 (' + lab4 + ') is on the wedge, pointing toward you. Trace anyway and flip.';
      else html = 'Group 4 (' + lab4 + ') is in the plane of the page. Swap it with the dash group, then trace.';
      $('swap').hidden = !(st.pose === 'plane' && !st.swaps);
    } else {
      if (!view) return;
      const o = R().orientation(st.coords, c.atom, order()[3], view.toViewer());
      four = o.four;
      html = four === 'back' ? 'Group 4 (' + lab4 + ') points away from you. Go to step 5.'
        : four === 'front' ? 'Group 4 (' + lab4 + ') points toward you. Tap Put 4 in back, or trace anyway and flip.'
          : 'Group 4 (' + lab4 + ') points sideways. Turn the model or tap Put 4 in back.';
    }
    const changed = four !== st.four;
    st.four = four;
    if (force || changed || $('orient-state').dataset.v !== html) { $('orient-state').innerHTML = html; $('orient-state').dataset.v = html; }
    if (st.stage === 'orient' || st.stage === 'answer') {
      const ok = traceAllowed();
      setDone('c-orient', ok);
      setCard('c-trace', ok);
      if (!ok && (st.trace.length || st.turn)) { clearTrace(); say('trace-msg', 'Group 4 is sideways now, so the trace was cleared.'); }
      if (ok && st.stage === 'orient' && !st.trace.length) $('trace-q').textContent = traceQ(four);
    }
    if (st.stage !== 'done') $('r-four').textContent = fourText(four) + (st.swaps ? ' (after a swap)' : '');
  }
  function traceQ(four) {
    const where = st.view === 'paper' ? 'on the drawing' : 'in the model';
    return 'Tap group 1, then 2, then 3, ' + where + ' or on the chips below' + (four === 'front' ? ' (4 points toward you, so flip at the end).' : '.');
  }
  function traceAllowed() {
    if (st.view === 'paper') return st.pose !== 'plane' || st.swaps > 0;
    return st.four === 'back' || st.four === 'front';
  }
  function putBack(front) {
    const c = cur();
    if (!c || !view || st.ranked.length < 4) return Promise.resolve();
    if (st.view !== '3d') setView('3d');
    if (st.stage === 'answer') clearTrace();
    if (st.tb) { st.tb = null; updateStyle(); }
    showViewerOnPhone();
    return view.orient({ center: c.atom, away: order()[3], up: order()[0], front: !!front }).then(() => { updateOrientation(true); });
  }
  function showViewerOnPhone() {
    if (narrow() && st.ready) document.querySelector('.viewer').scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' });
  }
  function setPose(pose, variant) {
    if (['dash', 'wedge', 'plane'].indexOf(pose) < 0) return;
    st.pose = pose; if (variant != null) st.variant = variant | 0;
    st.swaps = 0;
    ['dash', 'wedge', 'plane'].forEach(p => $('pose-' + p).setAttribute('aria-pressed', String(p === pose)));
    if (st.stage === 'answer') clearTrace(); else if (st.stage === 'orient') clearTrace();
    computeSlots(); updateOrientation(true); drawPaper();
  }
  function paperNew() {
    st.newCount++;
    st.variant++;
    if (st.newCount % 3 === 0) {
      const cyc = { dash: 'wedge', wedge: 'plane', plane: 'dash' };
      setPose(cyc[st.pose]);
    } else setPose(st.pose);
  }
  function doSwap() {
    if (st.view !== 'paper' || st.pose !== 'plane' || st.swaps) return;
    st.swaps = 1;
    const cv = $('paper'); cv.classList.remove('fade'); void cv.offsetWidth; if (!reduced()) cv.classList.add('fade');
    computeSlots(); updateOrientation(true); drawPaper();
  }

  /* ------------------------------------------------------------------ */
  /* step 5: trace                                                        */
  /* ------------------------------------------------------------------ */
  function clearTrace() {
    st.trace = []; st.turn = null; st.turnKnown = false; st.arrowT = 0; st.traceSnap = null; st.arrowHidden = false;
    $('turn').hidden = true; $('turn-row').hidden = true; $('answer-row').hidden = true; $('trace-tools').hidden = false; $('trace-chips').hidden = false;
    resetPickRows();
    if (st.stage === 'answer') st.stage = 'orient';
    syncTraceChips();
    drawPaper();
  }
  function resetPickRows() {
    ['ans-R', 'ans-S', 'turn-cw', 'turn-ccw'].forEach(id => { $(id).classList.remove('right', 'wrongans'); $(id).disabled = false; });
  }
  // step 5 chips: the ranked groups 1-3, a keyboard / phone way to trace
  function renderTraceChips() {
    const box = $('trace-chips'), c = cur();
    box.innerHTML = '';
    if (!c || st.ranked.length < 4) return;
    order().slice(0, 3).forEach((a, i) => {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'chip'; b.dataset.atom = a;
      b.innerHTML = '<span class="badge">' + (i + 1) + '</span><span class="gl">' + gHTML(a) + '</span>';
      b.setAttribute('aria-label', 'Group ' + (i + 1) + ', ' + gname(a));
      b.addEventListener('click', () => tapGroup(a));
      box.appendChild(b);
    });
    syncTraceChips();
  }
  function syncTraceChips() {
    document.querySelectorAll('#trace-chips .chip').forEach(b => b.setAttribute('aria-pressed', String(st.trace.indexOf(+b.dataset.atom) >= 0)));
  }
  function traceTap(X) {
    const ord = order();
    const k = st.trace.length;
    if (k >= 3) return;
    if (st.tb) st.tb = null;
    if (X === ord[k]) {
      st.trace.push(X);
      if (st.trace.length < 3) say('trace-msg', st.trace.length === 1 ? 'Good, that is 1. Now 2.' : 'Good, 2. Now 3.', 'ok');
      else completeTrace();
      syncTraceChips(); updateStyle(); drawPaper();
      return;
    }
    const n = ord.indexOf(X) + 1;
    if (n === 4) say('trace-msg', "That's group 4. Leave it out: trace only 1 → 2 → 3.", 'err');
    else say('trace-msg', k === 0 ? "That's group " + n + '. Start at 1.' : "That's group " + n + '. Next is ' + (k + 1) + '.', 'err');
    flash([X], 'wrong', 900);
    st.paperWrong = X; drawPaper();
    setTimeout(() => { if (st.paperWrong === X) { st.paperWrong = null; drawPaper(); } }, 900);
  }
  function screenPoints() {
    const c = cur(), ord = order();
    if (st.view === 'paper') {
      const g = PA().geometry(paperData(), paperCtx());
      return { c: g.center, p: [0, 1, 2].map(i => [g.groups[i].x, g.groups[i].y]) };
    }
    const pc = view.project(c.atom);
    return { c: [pc.x, pc.y], p: [0, 1, 2].map(i => { const q = view.project(ord[i]); return [q.x, q.y]; }) };
  }
  function completeTrace() {
    const sp = screenPoints();
    st.turn = R().turn2D(sp.c, sp.p[0], sp.p[1], sp.p[2]);
    if (!st.turn) { st.trace = []; syncTraceChips(); say('trace-msg', 'The three groups line up from here. Turn the model a little and trace again.', 'err'); return; }
    st.stage = 'answer';
    st.turnKnown = false;
    st.arrowT = reduced() ? 1 : 0;
    st.arrowStart = performance.now();
    if (view && st.view === '3d') {
      const f = view.screenFrame();
      st.traceSnap = { v: view.toViewer(), up: f.up };
    }
    const expect = R().labelFromTrace(st.turn, st.four, st.swaps);
    if (expect !== cur().label && !st.warned) { st.warned = true; console.warn('R/S: trace label', expect, 'differs from the engine label', cur().label); }
    // first the student reads the direction off the arrow, then picks R or S
    resetPickRows();
    $('turn').hidden = true;
    $('turn-row').hidden = false; $('answer-row').hidden = true; $('trace-tools').hidden = true; $('trace-chips').hidden = true;
    say('trace-msg', '');
    $('trace-q').textContent = 'The arrow is drawn. Follow it from 1 to 2 to 3.';
    drawPaper();
    reveal('turn-row');
  }
  function flipReason() { return st.swaps ? 'I swapped two groups' : st.four === 'front' ? '4 points toward you' : ''; }
  function turnPick(dir, shown) {
    if (st.stage !== 'answer' || st.turnKnown) return;
    const tw = turnWord(st.turn);
    if (!shown && dir !== st.turn) {
      $('turn-' + dir).classList.add('wrongans'); $('turn-' + dir).disabled = true;
      say('trace-msg', 'Not quite. Start at 1 and follow the arrow through 2 to 3. Is that the way a clock\'s hands move?', 'err');
      return;
    }
    st.turnKnown = true;
    $('turn-' + st.turn).classList.add('right');
    $('turn-cw').disabled = true; $('turn-ccw').disabled = true;
    $('turn').textContent = cap(tw); $('turn').hidden = false;
    $('turn-row').hidden = true; $('answer-row').hidden = false;
    const why = flipReason();
    say('trace-msg', (shown ? 'It turns ' + tw + '.' : 'Yes, ' + tw + '.') + (why ? ' Remember: ' + why + '.' : '') + ' So is it R or S?', shown ? '' : 'ok');
    $('trace-q').textContent = 'Now name it.';
    drawPaper();
    reveal('answer-row');
  }
  function rightSentence() {
    const lab = cur().label, tw = turnWord(st.turn);
    if (st.swaps) return 'Right. The arrow goes ' + tw + ', but I swapped two groups, so flip it: ' + lab + '.';
    if (st.four === 'front') return 'Right. The arrow goes ' + tw + ', but 4 points toward you, so flip it: ' + lab + '.';
    return 'Right. ' + cap(tw) + ' with 4 in back is ' + lab + '.';
  }
  function answer(L, shown) {
    if (st.stage !== 'answer') return;
    if (!st.turnKnown) turnPick(st.turn, true);   // "Show me" or the test handle skips the direction question
    const c = cur(), lab = c.label;
    const why = flipReason();
    if (!shown && L !== lab) {
      // one more try: say it's wrong, give a hint, keep the right letter hidden
      if (L === 'R' || L === 'S') { $('ans-' + L).classList.add('wrongans'); $('ans-' + L).disabled = true; }
      st.ansMissed = true;
      say('trace-msg', why
        ? 'Not quite. ' + cap(why) + ', so the letter is the opposite of what the arrow says. Try again.'
        : 'Not quite. With 4 in back, read the arrow as it is: clockwise is R, counterclockwise is S. Try again.', 'err');
      return;
    }
    const msg = shown ? "Here's how it goes: " + rightSentence().replace(/^Right\. /, '') : rightSentence();
    say('trace-msg', msg, 'ok');
    $('ans-' + lab).classList.add('right');
    $('ans-R').disabled = true; $('ans-S').disabled = true;
    st.answers[c.atom] = lab;
    st.stage = 'done';
    $('r-label').textContent = lab;
    $('r-four').textContent = st.swaps ? '4 swapped onto the dash, so flip' : st.four === 'front' ? '4 toward you, so flip' : '4 away from you';
    setDone('c-trace', true);
    if (view) view.setPickable([]);
    const left = st.centers.length - Object.keys(st.answers).length;
    $('r-summary').textContent = st.centers.length > 1 ? (st.centers.length - left) + ' of ' + st.centers.length + ' stereocenters done' : '1 stereocenter';
    if (left > 0) { $('next-center').hidden = false; reveal('next-center'); }
    else finish();
    updateStyle(); drawPaper();
  }
  function traceShow() {
    if (st.stage !== 'orient') return Promise.resolve();
    const go = () => { st.trace = []; [0, 1, 2].forEach(i => traceTap(order()[i])); };
    if (!traceAllowed()) {
      if (st.view === 'paper') { if (st.pose === 'plane') doSwap(); else setPose('dash'); go(); return Promise.resolve(); }
      return putBack(false).then(go);
    }
    go();
    return Promise.resolve();
  }
  function nextCenter() {
    const k = st.centers.findIndex(c => !st.answers[c.atom]);
    if (k < 0) return;
    startCenter(k);
    if (narrow()) $('c-rank').scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' });
  }

  /* ------------------------------------------------------------------ */
  /* step 6: name                                                         */
  /* ------------------------------------------------------------------ */
  function finish() {
    st.stage = 'finished';
    const labels = {};
    st.centers.forEach(c => { labels[c.atom] = c.label; });
    let opts = {}, numberingNote = '';
    if (st.preset) opts = { locants: st.preset.locants, base: st.preset.base };
    else {
      const loc = {}; let all = true;
      st.centers.forEach(c => { const m = /^C(\d+)$/.exec(centerName(c.atom)); if (m) loc[c.atom] = +m[1]; else all = false; });
      let base = cleanName(st.mol.name || (st.src && st.src.name));
      if (!base) {
        try { const mc = C().mainChain(st.src || st.mol); if (mc && mc.supported && mc.name) base = cleanName(mc.name); } catch (e) { /* no name */ }
      }
      opts = all ? { locants: loc } : {};
      if (base) opts.base = base;
      numberingNote = all ? 'C numbers come from the main chain I found.' : '';
    }
    let d = R().describe(st.mol, labels, opts);
    if (!st.preset && d.meso && opts.locants) {
      // meso: numbering from either end is allowed, and IUPAC then gives R the lower locant
      const inv = {}; Object.keys(labels).forEach(a => { inv[a] = flip(labels[a]); });
      const d2 = R().describe(st.mol, inv, opts);
      const letters = x => x.descriptor.replace(/[^RS]/g, '');
      if (letters(d2) < letters(d)) { d = d2; numberingNote = 'This is meso, so the chain can be numbered from either end. IUPAC numbers it so R gets the lower number, which is the other end from the C numbers above.'; }
    }
    $('summary-text').textContent = d.text || d.descriptor;
    $('summary-text').classList.remove('wait');
    const notes = [];
    const matches = st.preset && Object.keys(st.preset.expect).every(a => labels[a] === st.preset.expect[a]);
    if (d.meso && !(st.preset && /^Meso:/.test(st.preset.note))) notes.push('Meso: it has stereocenters, but its mirror image is the same molecule.');
    if (st.preset && st.preset.note && (matches || st.preset.meso || /L\/D and R\/S/.test(st.preset.note))) notes.push(st.preset.note);
    if (st.mirrored && st.preset && !d.meso) notes.push('This is the mirror image of the preset, so every letter flipped.');
    if (numberingNote) notes.push(numberingNote);
    $('summary-note').textContent = notes.join(' ');
    $('r-summary').textContent = (d.text || d.descriptor) + (d.meso ? ' · meso' : '');
    setCard('c-name', true); setDone('c-name', true);
    $('next-center').hidden = true;
    reveal('c-name');
    syncDock();
    // phones: the answer shows in the result card under the model; then bring up the name
    if (narrow() && st.ready) setTimeout(() => { if (st.stage === 'finished') $('c-name').scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' }); }, reduced() ? 0 : 1200);
  }

  /* ------------------------------------------------------------------ */
  /* mirror / restart                                                     */
  /* ------------------------------------------------------------------ */
  function mirror() {
    if (!st.mol) return;
    st.coords = R().mirror(st.coords);
    st.mirrored = !st.mirrored;
    if (view) view.setCoords(st.coords);
    const found = new Set(st.found), wasFind = st.stage === 'find';
    computeCenters();
    st.answers = {};
    resetCenterSteps();
    ['c-rank', 'c-orient', 'c-trace', 'c-name'].forEach(id => { setCard(id, false); setDone(id, false); });
    $('summary-text').textContent = 'Finish each stereocenter first.'; $('summary-text').classList.add('wait'); $('summary-note').textContent = '';
    st.found = found; syncAtomChips();
    if (wasFind || st.stage === 'none' || !st.centers.length) {
      st.stage = 'find'; setCard('c-find', true); setDone('c-find', false); $('result').hidden = true;
    } else {
      st.stage = 'find'; startCenter(0);
    }
    say(msgId(), 'Mirror image on screen. Every R becomes S and every S becomes R.', 'ok');
    updateStyle(); drawPaper();
  }
  function restart() {
    if (!st.mol) return;
    if (view) view.reset();
    resetAll();
    if (narrow()) $('c-find').scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' });
  }

  /* ------------------------------------------------------------------ */
  /* view switching, styles, overlay, paper                               */
  /* ------------------------------------------------------------------ */
  function setView(v) {
    if (v === 'paper' && !cur()) return;
    st.view = v;
    $('view-3d').setAttribute('aria-pressed', String(v === '3d'));
    $('view-paper').setAttribute('aria-pressed', String(v === 'paper'));
    $('paper').hidden = v !== 'paper';
    document.querySelector('.viewer').classList.toggle('paper', v === 'paper');
    $('orient-3d').hidden = v !== '3d';
    $('orient-paper').hidden = v !== 'paper';
    if (st.stage === 'answer' || (st.stage === 'orient' && st.trace.length)) { clearTrace(); say('trace-msg', 'New view, so trace again.'); }
    else if (st.stage === 'done' || st.stage === 'finished') { st.turn = null; st.trace = []; }
    if (v === 'paper') { computeSlots(); ['dash', 'wedge', 'plane'].forEach(p => $('pose-' + p).setAttribute('aria-pressed', String(p === st.pose))); }
    updateOrientation(true);
    drawPaper();
  }

  function updateStyle() {
    if (!view || !st.mol) return;
    const map = {};
    const c = cur();
    if (st.stage === 'find' || !c) {
      st.found.forEach(a => { map[a] = 'found'; });
      if (st.findPair) { st.findPair.a.forEach(a => { map[a] = 'tbA'; }); st.findPair.b.forEach(a => { map[a] = 'tbB'; }); }
    } else {
      st.centers.forEach(x => { if (st.answers[x.atom]) map[x.atom] = 'found'; });
      c.ranking.ligands.forEach(l => { map[l.atom] = 'ligand'; });
      if (st.tb && st.stage === 'rank') { st.tb.a.forEach(a => { map[a] = 'tbA'; }); st.tb.b.forEach(a => { map[a] = 'tbB'; }); }
      if (st.tb && st.stage !== 'rank') { st.tb.a.forEach(a => { map[a] = 'tbA'; }); st.tb.b.forEach(a => { map[a] = 'tbB'; }); }
      map[c.atom] = 'center';
    }
    Object.keys(st.flash).forEach(a => { map[a] = st.flash[a].style; });
    view.setStyle(map);
    const tb = document.querySelector('#legend .lg-tb');
    if (tb) tb.hidden = !Object.keys(map).some(a => map[a] === 'tbA' || map[a] === 'tbB');
    syncDock();
  }

  function sizeCanvas(cv) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = cv.clientWidth || 1, h = cv.clientHeight || 1;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }
  const angleOf = (c, p) => Math.atan2(-(p[1] - c[1]), p[0] - c[0]) * 180 / Math.PI;

  function drawOverlay() {
    const cv = $('overlay');
    const { ctx, w, h } = sizeCanvas(cv);
    ctx.clearRect(0, 0, w, h);
    if (!view || !st.mol || st.view !== '3d') return;
    const c = cur();
    let arcTop = null;
    if (!c || st.stage === 'find' || (st.stage === 'finished' && !st.turn)) { drawTags(ctx, null); return; }
    const ord = order();
    // trace path
    if (st.trace.length >= 2) {
      const pts = st.trace.map(a => view.project(a));
      ctx.lineCap = 'round'; ctx.lineJoin = 'round';
      [[6, '#0f0b13'], [3, '#ffb300']].forEach(([lw, col]) => {
        ctx.setLineDash(lw === 3 ? [8, 6] : []);
        ctx.strokeStyle = col; ctx.lineWidth = lw; ctx.beginPath();
        pts.forEach((p, i) => (i ? ctx.lineTo(p.x, p.y) : ctx.moveTo(p.x, p.y))); ctx.stroke();
      });
      ctx.setLineDash([]);
    }
    // arc radius around the center, the widest free gap for badge 4 (when it hides behind the center) and the labels
    const pc0 = view.project(c.atom), cc = [pc0.x, pc0.y];
    const norm360 = x => ((x % 360) + 360) % 360;
    const widestGap = angs => {
      const a = angs.map(norm360).sort((x, y) => x - y);
      if (!a.length) return 90;
      let best = -1, at = 90;
      a.forEach((x, i) => { const nx = i + 1 < a.length ? a[i + 1] : a[0] + 360; if (nx - x > best) { best = nx - x; at = x + (nx - x) / 2; } });
      return norm360(at);
    };
    let ar = 0, b4 = null, labelAt = 90;
    if (st.ranked.length === 4) {
      const pts = [0, 1, 2].map(i => { const q = view.project(ord[i]); return [q.x, q.y]; });
      const rr = pts.reduce((s, p) => s + Math.hypot(p[0] - cc[0], p[1] - cc[1]), 0) / 3;
      ar = Math.max(26, rr * 0.6);
      const angs = pts.filter(p => Math.hypot(p[0] - cc[0], p[1] - cc[1]) >= 6).map(p => angleOf(cc, p));
      const q4 = view.project(ord[3]), L4 = Math.hypot(q4.x - cc[0], q4.y - cc[1]);
      if (L4 < ar * 0.6) b4 = widestGap(angs);                     // 4 points (nearly) straight back
      else angs.push(angleOf(cc, [q4.x, q4.y]));
      labelAt = widestGap(b4 == null ? angs : angs.concat([b4]));
      // arrow
      if (st.turn && (st.stage === 'answer' || st.stage === 'done' || st.stage === 'finished') && !st.arrowHidden) {
        const answered = !!st.answers[c.atom];
        // the direction word waits until the student has said it
        const rad = labelAt * Math.PI / 180;
        PA().drawArrow(ctx, { cx: cc[0], cy: cc[1], r: ar, from: angleOf(cc, pts[0]), via: angleOf(cc, pts[1]), to: angleOf(cc, pts[2]), dir: st.turn, t: st.arrowT,
          textAt: labelAt, textR: ar + 14 + Math.abs(Math.cos(rad)) * (st.turn === 'cw' ? 33 : 58) + Math.abs(Math.sin(rad)) * 9,
          text: !answered && !!st.turnKnown });
        if (answered) arcTop = [cc[0] + Math.cos(rad) * (ar + 26), cc[1] - Math.sin(rad) * (ar + 26)];
      }
    }
    // badges, pushed out from the center along each bond; a hidden group 4 gets its badge in the free gap, clear of the arc
    st.ranked.forEach((a, k) => {
      const p = view.project(a), r = view.projRadius(a);
      if (k === 3 && b4 != null) {
        const rad = b4 * Math.PI / 180;
        PA().drawBadge(ctx, cc[0] + Math.cos(rad) * (ar + 22), cc[1] - Math.sin(rad) * (ar + 22), '4');
        return;
      }
      let dx = p.x - pc0.x, dy = p.y - pc0.y;
      const L = Math.hypot(dx, dy);
      if (L < 8) { dx = 0.7; dy = 0.7; } else { dx /= L; dy /= L; }
      PA().drawBadge(ctx, p.x + dx * (r + 13), p.y + dy * (r + 13), String(k + 1));
    });
    drawTags(ctx, arcTop);
  }
  function drawTags(ctx, arcTop) {
    // R/S tags on answered centers (the current one sits above its arrow when the arrow is shown)
    Object.keys(st.answers).forEach(a => {
      a = +a;
      const p = view.project(a), r = view.projRadius(a);
      let x = p.x + r + 16, y = p.y - r - 14;
      if (arcTop && cur() && a === cur().atom) { x = arcTop[0]; y = arcTop[1]; }
      ctx.beginPath(); ctx.arc(x, y, 17, 0, Math.PI * 2); ctx.fillStyle = 'rgba(23,19,29,.92)'; ctx.fill();
      ctx.lineWidth = 2.5; ctx.strokeStyle = '#ffb300'; ctx.stroke();
      ctx.font = '400 21px Shrikhand, Georgia, serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillStyle = '#ffb300'; ctx.fillText(st.answers[a], x, y + 1);
    });
  }

  let paperCtxCache = null;
  function paperCtx() { if (!paperCtxCache) paperCtxCache = $('paper').getContext('2d'); return paperCtxCache; }
  function paperLabel(a) {
    return gname(a).replace('(toward ', '(to ').replace(' (other side)', ' (other)');
  }
  function paperData() {
    const cv = $('paper'), c = cur();
    const w = cv.clientWidth || 300, h = cv.clientHeight || 300;
    const ord = c.ranking.order;
    if (!st.slots) computeSlots();
    const chem = C();
    const groups = ord.map((a, i) => {
      const el = st.mol.atoms[a].el;
      const k = st.ranked.indexOf(a);
      let state = 'idle';
      if (st.paperWrong === a) state = 'wrong';
      else if (st.trace.indexOf(a) >= 0) state = 'ok';
      return { label: paperLabel(a), color: (chem.EL[el] || chem.EL.C).color2d === '#ffb300' ? '#fbf3e2' : (chem.EL[el] || chem.EL.C).color2d,
        slot: st.slots[i], badge: k >= 0 ? String(k + 1) : null, state };
    });
    let arrow = null;
    if (st.turn && st.view === 'paper') {
      const deg = i => PA().SLOTS[st.slots[i]].deg;
      arrow = { from: deg(0), via: deg(1), to: deg(2), dir: st.turn, t: st.arrowT, text: !!st.turnKnown };
    }
    const cap4 = st.ranked.length < 4 ? (st.stage === 'rank' ? 'Tap the groups from 1 to 4.' : '')
      : st.swaps ? 'Group 4 swapped onto the dash (mirror image)' : st.pose === 'dash' ? 'Group 4 on the dash, away from you'
        : st.pose === 'wedge' ? 'Group 4 on the wedge, toward you' : 'Group 4 in the plane of the page';
    const shiftX = narrow() ? 0 : Math.round(Math.min(110, w * 0.1));
    return { w, h, shiftX, centerLabel: centerName(c.atom), groups, arrow, caption: cap4, tag: st.answers[c.atom] || null };
  }
  function drawPaper() {
    const cv = $('paper');
    if (st.view !== 'paper' || !cur()) return;
    const { ctx } = sizeCanvas(cv);
    const data = paperData();
    PA().draw(ctx, data);
    const where = { left: 'on the left', right: 'on the right', wedge: 'on a wedge', dash: 'on a dash' };
    cv.setAttribute('aria-label', 'Stereocenter ' + data.centerLabel + '. ' + data.groups.map(g => g.label + ' ' + where[g.slot]).join(', ') + '.');
  }

  function renderLegend() {
    const chem = C(), seen = [];
    st.mol.atoms.forEach(a => { if (seen.indexOf(a.el) < 0) seen.push(a.el); });
    const ord = ['C', 'H', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'B', 'P'];
    seen.sort((a, b) => ord.indexOf(a) - ord.indexOf(b));
    const tip = coarse() ? 'drag to rotate · tap an atom' : 'drag to rotate · click an atom';
    $('legend').innerHTML = seen.map(el => {
      const c = '#' + (chem.EL[el] || chem.EL.C).color3d.toString(16).padStart(6, '0');
      return '<span><i style="background:' + c + ';box-shadow:0 0 0 1px rgba(251,243,226,.45)"></i>' + el + '</span>';
    }).join('') + '<span><i class="lg-ring"></i>stereocenter</span>' +
      '<span class="lg-tb" hidden><i style="background:#14b8a6"></i><i style="background:#ff5c93"></i>groups being compared</span>' +
      '<span class="tip">' + tip + '</span>';
  }

  /* ------------------------------------------------------------------ */
  /* frame loop hooks                                                     */
  /* ------------------------------------------------------------------ */
  function onViewChange() { st.orientDirty = true; }
  function onFrame() {
    if (st.turn && st.arrowT < 1) {
      st.arrowT = Math.min(1, (performance.now() - (st.arrowStart || 0)) / 600);
      if (st.view === 'paper') drawPaper();
    }
    if (st.orientDirty && view && !view.busy()) {
      st.orientDirty = false;
      if (st.view === '3d') {
        // a real camera change after the arrow is drawn clears the trace
        if (st.traceSnap && st.turn) {
          const f = view.screenFrame(), v = view.toViewer();
          const d = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
          const moved = d(v, st.traceSnap.v) < 0.9976 || d(f.up, st.traceSnap.up) < 0.9976;   // about 4 degrees
          if (moved) {
            if (st.stage === 'answer') { clearTrace(); say('trace-msg', 'You moved the model, so trace again.'); }
            else st.arrowHidden = true;
            st.traceSnap = null;
          }
        }
        updateOrientation(false);
      }
    }
    drawOverlay();
  }

  /* ------------------------------------------------------------------ */
  /* input modes, sketch, presets, hash                                   */
  /* ------------------------------------------------------------------ */
  function setInputMode(m) {
    ['type', 'draw', 'presets'].forEach(k => {
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
    const sel = $('preset'), P = PR();
    sel.innerHTML = '<option value="" disabled selected>Choose a molecule…</option>' + P.groups.map(g =>
      '<optgroup label="' + esc(g) + '">' + P.list.filter(p => p.group === g).map(p =>
        '<option value="' + p.id + '">' + esc(p.name) + '</option>').join('') + '</optgroup>').join('');
  }
  function fromHash() {
    let h = '';
    try { h = decodeURIComponent((location.hash || '').slice(1)); } catch (e) { h = (location.hash || '').slice(1); }
    if (!h || '#' + encodeURIComponent(h) === lastHash || '#' + h === lastHash) return false;
    const p = PR().byId(h);
    if (p) return loadPreset(p, { hash: false, msgId: msgId() });
    const text = h.startsWith('smiles=') ? h.slice(7) : h;
    $('mol-input').value = text;
    return loadText(text, { hash: false, msgId: msgId() });
  }

  function wire() {
    $('mode-type').addEventListener('click', () => setInputMode('type'));
    $('mode-draw').addEventListener('click', () => setInputMode('draw'));
    $('mode-presets').addEventListener('click', () => setInputMode('presets'));
    $('mol-go').addEventListener('click', () => loadText($('mol-input').value));
    $('mol-input').addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); loadText($('mol-input').value); } });
    document.querySelectorAll('.hint .ex').forEach(b => b.addEventListener('click', () => { $('mol-input').value = b.textContent; loadText(b.textContent); }));
    $('preset').addEventListener('change', e => { const p = PR().byId(e.target.value); if (p) loadPreset(p, { msgId: 'preset-msg' }); });

    $('sk-free').addEventListener('click', () => setSketchMode('free'));
    $('sk-chain').addEventListener('click', () => setSketchMode('chain'));
    $('sk-finish').addEventListener('click', () => { if (sketch) sketch.setMode('chain'); });
    $('sk-undo').addEventListener('click', () => sketch && sketch.undo());
    $('sk-clear').addEventListener('click', () => sketch && sketch.clear());
    $('sk-build').addEventListener('click', () => {
      if (!sketch) return;
      const r = sketch.getMolecule();
      if (!r.ok) { say('sk-msg', r.error || 'Check the drawing.', 'err'); return; }
      if (loadMol(r.mol, { msgId: 'sk-msg' }) && narrow()) document.querySelector('.viewer').scrollIntoView({ behavior: reduced() ? 'auto' : 'smooth', block: 'start' });
    });
    $('mirror').addEventListener('click', mirror);

    $('find-done').addEventListener('click', findDone);
    $('find-none').addEventListener('click', findNone);
    $('find-show').addEventListener('click', findShow);
    $('rank-undo').addEventListener('click', rankUndo);
    $('rank-show').addEventListener('click', rankShow);
    $('put-back').addEventListener('click', () => putBack(false));
    $('put-front').addEventListener('click', () => putBack(true));
    ['dash', 'wedge', 'plane'].forEach(p => $('pose-' + p).addEventListener('click', () => { setPose(p); showViewerOnPhone(); }));
    $('paper-new').addEventListener('click', () => { paperNew(); showViewerOnPhone(); });
    $('swap').addEventListener('click', () => { doSwap(); showViewerOnPhone(); });
    $('view-paper').addEventListener('click', showViewerOnPhone);
    $('ans-R').addEventListener('click', () => answer('R'));
    $('ans-S').addEventListener('click', () => answer('S'));
    $('ans-show').addEventListener('click', () => answer(null, true));
    $('turn-cw').addEventListener('click', () => turnPick('cw'));
    $('turn-ccw').addEventListener('click', () => turnPick('ccw'));
    $('turn-show').addEventListener('click', () => turnPick(null, true));
    $('trace-go').addEventListener('click', showViewerOnPhone);
    if (window.matchMedia) { const mq = window.matchMedia('(max-width: 860px)'); const f = () => syncDock(); if (mq.addEventListener) mq.addEventListener('change', f); else if (mq.addListener) mq.addListener(f); }
    $('trace-show').addEventListener('click', traceShow);
    $('next-center').addEventListener('click', nextCenter);
    $('restart').addEventListener('click', restart);
    $('view-3d').addEventListener('click', () => setView('3d'));
    $('view-paper').addEventListener('click', () => setView('paper'));

    // paper taps
    const pc = $('paper');
    let down = null;
    pc.addEventListener('pointerdown', e => { down = { x: e.clientX, y: e.clientY }; });
    pc.addEventListener('pointerup', e => {
      if (!down || Math.hypot(e.clientX - down.x, e.clientY - down.y) > 8) { down = null; return; }
      down = null;
      const r = pc.getBoundingClientRect();
      if (!cur()) return;
      const i = PA().hit(paperData(), e.clientX - r.left, e.clientY - r.top, paperCtx());
      if (i == null) return;
      tapGroup(order()[i]);
    });
    pc.addEventListener('pointermove', e => {
      if (e.pointerType !== 'mouse' || !cur()) return;
      const r = pc.getBoundingClientRect();
      pc.style.cursor = PA().hit(paperData(), e.clientX - r.left, e.clientY - r.top, paperCtx()) != null ? 'pointer' : '';
    });
    new ResizeObserver(() => drawPaper()).observe(pc);
    window.addEventListener('hashchange', fromHash);
  }

  /* ------------------------------------------------------------------ */
  /* start                                                                */
  /* ------------------------------------------------------------------ */
  function start() {
    if (!window.NNChem || !window.NNGeom || !window.RSCip || !window.RSPresets || !window.RSPaper || !window.RSView3D) {
      say('preset-msg', 'Part of this tool did not load. Reload the page to try again.', 'err');
      return;
    }
    const host = $('canvas');
    const forceFlat = /[?&]flat=1/.test(location.search);
    view = RSView3D.create(host, { onTap: i => tapAtom(i), onChange: onViewChange, onFrame, forceFlat });
    if (!view) {
      const d = document.createElement('div'); d.className = 'nogl';
      d.textContent = "Your browser can't draw the model here. The Paper view and the steps still work.";
      host.appendChild(d);
    }
    fillPresets();
    wire();
    if (!fromHash()) loadPreset(PR().byId('r-2-butanol'), { hash: false, msgId: 'preset-msg' });
    st.ready = true;
    if (!view) setInterval(onFrame, 100);
  }

  window.RSApp = {
    load: text => loadText(String(text)),
    loadPreset: id => loadPreset(PR().byId(id), {}),
    get state() { return st; },
    get view() { return view; },
    tapAtom: i => tapAtom(i),
    tapGroup: a => tapGroup(a),
    setView: v => setView(v),
    setPose: (p, v) => setPose(p, v),
    putBack: () => putBack(false),
    putFront: () => putBack(true),
    answer: L => answer(L),
    turn: d => turnPick(d),
    paperData: () => (cur() ? paperData() : null),
    findDone, findNone, findShow, rankShow, rankUndo, traceShow, next: nextCenter, mirror, swap: doSwap
  };
  window.addEventListener('load', start);
})();

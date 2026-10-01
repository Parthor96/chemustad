/* rs/presets.js: window.RSPresets, teaching molecules for the R/S Assigner and a name lookup.
 * Pure, no DOM, node-requirable. Contract: rs/SPEC.md section 4.
 * Every SMILES label here is RDKit-verified (rs/test/rdkit-cases.json). Atom numbers in `locants` and
 * `expect` are heavy-atom indices in the SMILES.
 */
(function (root) {
  'use strict';

  const NOTE_LD = 'L/D and R/S are separate systems. They often line up, but not always.';
  const NOTE_MESO = 'Meso: it has stereocenters, but its mirror image is the same molecule.';

  function P(id, name, aliases, smiles, base, locants, expect, meso, group, note) {
    return { id, name, aliases, smiles, base, locants, group, expect, meso: !!meso, note: note || '' };
  }

  const G1 = 'One stereocenter', GT = 'Tie-breaks', GD = 'Double & triple bonds', GR = 'Rings',
    G2 = 'Two stereocenters', GN = 'No stereocenter';

  const list = [
    P('r-2-butanol', '(R)-2-butanol', ['2-butanol', 'butan-2-ol', 'sec-butyl alcohol'], 'C[C@@H](O)CC', '2-butanol',
      { 1: 2 }, { 1: 'R' }, false, G1, ''),
    P('s-2-bromobutane', '(S)-2-bromobutane', ['2-bromobutane'], 'C[C@H](Br)CC', '2-bromobutane',
      { 1: 2 }, { 1: 'S' }, false, G1, ''),
    P('r-2-chlorobutane', '(R)-2-chlorobutane', ['2-chlorobutane'], 'C[C@@H](Cl)CC', '2-chlorobutane',
      { 1: 2 }, { 1: 'R' }, false, G1, ''),
    P('r-bcf', '(R)-bromochlorofluoromethane', ['bromochlorofluoromethane', 'CHBrClF'], 'F[C@H](Cl)Br', 'bromochlorofluoromethane',
      null, { 1: 'R' }, false, G1, 'Every group starts with a different atom, so no tie-breaks are needed.'),
    P('s-alanine', '(S)-alanine', ['L-alanine', 'alanine'], 'C[C@H](N)C(=O)O', 'alanine',
      { 1: 2 }, { 1: 'S' }, false, G1, NOTE_LD),
    P('s-lactic', '(S)-lactic acid', ['L-lactic acid', 'lactic acid'], 'C[C@H](O)C(=O)O', 'lactic acid',
      { 1: 2 }, { 1: 'S' }, false, G1, ''),
    P('r-glyceraldehyde', '(R)-glyceraldehyde', ['D-glyceraldehyde', 'glyceraldehyde'], 'O=C[C@H](O)CO', 'glyceraldehyde',
      { 2: 2 }, { 2: 'R' }, false, G1, NOTE_LD),

    P('s-3-methylhexane', '(S)-3-methylhexane', ['3-methylhexane'], 'CC[C@@H](CCC)C', '3-methylhexane',
      { 2: 3 }, { 2: 'S' }, false, GT, ''),
    P('s-trimethylbutanol', '(S)-2,3,3-trimethyl-1-butanol', ['2,3,3-trimethylbutan-1-ol', '2,3,3-trimethyl-1-butanol'],
      'OC[C@@H](C)C(C)(C)C', '2,3,3-trimethyl-1-butanol', { 2: 2 }, { 2: 'S' }, false, GT,
      'CH2OH beats C(CH3)3: compare (O,H,H) with (C,C,C). O wins at the first position of the set.'),
    P('r-chlorodimethylbutane', '(R)-1-chloro-2,3-dimethylbutane', ['1-chloro-2,3-dimethylbutane'], 'ClC[C@H](C)C(C)C',
      '1-chloro-2,3-dimethylbutane', { 2: 2 }, { 2: 'R' }, false, GT, ''),
    P('r-cysteine', '(R)-cysteine', ['L-cysteine', 'cysteine'], 'N[C@@H](CS)C(=O)O', 'cysteine',
      { 1: 2 }, { 1: 'R' }, false, GT,
      'Most natural amino acids are S. Cysteine is R because sulfur (16) outranks the oxygens of COOH. L and D are a different naming system from R and S.'),

    P('r-butenol', '(R)-3-buten-2-ol', ['but-3-en-2-ol', '3-buten-2-ol'], 'C[C@@H](O)C=C', '3-buten-2-ol',
      { 1: 2 }, { 1: 'R' }, false, GD, ''),
    P('s-dimethylpentene', '(S)-3,4-dimethyl-1-pentene', ['3,4-dimethylpent-1-ene', '3,4-dimethyl-1-pentene'], 'C=C[C@@H](C)C(C)C',
      '3,4-dimethyl-1-pentene', { 2: 3 }, { 2: 'S' }, false, GD, 'The double bond counts its carbon twice.'),
    P('s-trimethylpentyne', '(S)-3,4,4-trimethyl-1-pentyne', ['3,4,4-trimethylpent-1-yne', '3,4,4-trimethyl-1-pentyne'],
      'C#C[C@@H](C)C(C)(C)C', '3,4,4-trimethyl-1-pentyne', { 2: 3 }, { 2: 'S' }, false, GD, 'The triple bond counts its carbon three times.'),

    P('r-limonene', '(R)-limonene', ['limonene', 'D-limonene'], 'CC1=CC[C@@H](CC1)C(=C)C', 'limonene',
      { 4: 4 }, { 4: 'R' }, false, GR, ''),
    P('r-methylcyclohexanone', '(R)-3-methylcyclohexanone', ['3-methylcyclohexanone'], 'C[C@@H]1CCCC(=O)C1', '3-methylcyclohexanone',
      { 1: 3 }, { 1: 'R' }, false, GR, ''),
    P('cis-12-dmch', 'cis-1,2-dimethylcyclohexane', ['cis-1,2-dimethylcyclohexane'], 'C[C@@H]1CCCC[C@@H]1C', '1,2-dimethylcyclohexane',
      { 1: 1, 6: 2 }, { 1: 'R', 6: 'S' }, true, GR, NOTE_MESO),
    P('trans-12-dmch', 'trans-1,2-dimethylcyclohexane', ['trans-1,2-dimethylcyclohexane'], 'C[C@@H]1CCCC[C@H]1C', '1,2-dimethylcyclohexane',
      { 1: 1, 6: 2 }, { 1: 'R', 6: 'R' }, false, GR, ''),
    P('cis-13-dmch', 'cis-1,3-dimethylcyclohexane', ['cis-1,3-dimethylcyclohexane'], 'C[C@@H]1CCC[C@H](C)C1', '1,3-dimethylcyclohexane',
      { 1: 1, 5: 3 }, { 1: 'R', 5: 'S' }, true, GR, NOTE_MESO),
    P('trans-2-mch', 'trans-2-methylcyclohexanol', ['trans-2-methylcyclohexanol'], 'C[C@@H]1CCCC[C@H]1O', '2-methylcyclohexanol',
      { 6: 1, 1: 2 }, { 6: 'R', 1: 'R' }, false, GR, ''),

    P('meso-dibromobutane', 'meso-2,3-dibromobutane', ['(2R,3S)-2,3-dibromobutane'], 'C[C@@H](Br)[C@@H](Br)C', '2,3-dibromobutane',
      { 1: 2, 3: 3 }, { 1: 'R', 3: 'S' }, true, G2, NOTE_MESO),
    P('rr-dibromobutane', '(2R,3R)-2,3-dibromobutane', [], 'C[C@@H](Br)[C@H](Br)C', '2,3-dibromobutane',
      { 1: 2, 3: 3 }, { 1: 'R', 3: 'R' }, false, G2, ''),
    P('methylhexanol', '(2R,3S)-3-methyl-2-hexanol', ['3-methyl-2-hexanol', '3-methylhexan-2-ol'], 'C[C@@H](O)[C@@H](C)CCC',
      '3-methyl-2-hexanol', { 1: 2, 3: 3 }, { 1: 'R', 3: 'S' }, false, G2, ''),
    P('threonine', '(2S,3R)-threonine', ['L-threonine', 'threonine'], 'C[C@@H](O)[C@H](N)C(=O)O', 'threonine',
      { 3: 2, 1: 3 }, { 3: 'S', 1: 'R' }, false, G2, ''),

    P('3-methylpentane', '3-methylpentane', [], 'CCC(C)CC', '3-methylpentane',
      null, {}, false, GN, 'C3 looks busy, but two of its groups are ethyl.'),
    P('cis-4-mch', 'cis-4-methylcyclohexanol', ['cis-4-methylcyclohexanol'], 'C[C@H]1CC[C@@H](O)CC1', '4-methylcyclohexanol',
      null, {}, false, GN, 'Cis and trans isomers exist here. They are named cis or trans, not R/S.')
  ];

  const groups = [G1, GT, GD, GR, G2, GN];

  function byId(id) { return list.find(p => p.id === id) || null; }

  const norm = t => String(t || '').toLowerCase().replace(/[−–—]/g, '-').replace(/\s+/g, ' ').trim();
  const DESC = /^\(?((?:\d?[rs])(?:,\s*\d?[rs])*)\)?-/i;

  function normDescriptor(d) {
    const parts = d.split(',').map(s => s.trim().toUpperCase()).filter(Boolean);
    return '(' + parts.join(',') + ')';
  }

  function lookup(text) {
    const t = norm(text);
    if (!t) return null;
    const names = p => [p.name, p.base].concat(p.aliases).map(norm);
    // exact match with the descriptor still attached (meso-..., (2R,3S)-..., cis-...)
    let hit = list.find(p => names(p).indexOf(t) >= 0);
    if (hit) return { preset: hit, want: null };
    const m = t.match(DESC);
    if (!m) return null;
    const rest = t.slice(m[0].length).trim();
    hit = list.find(p => names(p).indexOf(rest) >= 0);
    if (!hit) return null;
    return { preset: hit, want: normDescriptor(m[1]) };
  }

  const api = { list, groups, byId, lookup, _norm: norm, DESC };
  root.RSPresets = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);

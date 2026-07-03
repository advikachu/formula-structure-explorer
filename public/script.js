/* =================================================================
   CHEMICAL FORMULA -> STRUCTURE / HYBRIDIZATION / BOND ANGLE TOOL
   -----------------------------------------------------------------
   Pipeline:
     1. Parse the typed formula into element counts.
     2. Look it up in a built-in library of common compounds
        (each stored as a tiny SMILES string + friendly name).
        If multiple known compounds share that formula (isomers),
        ask the user to choose one.
     3. If not in the library, try a "single central atom" VSEPR
        heuristic (covers things like CO2, NH3, SF6, PCl5, BF3...).
     4. Build a molecule graph (atoms + bonds + implicit H),
        lay it out in 2D with a force-directed algorithm, draw it,
        and compute per-atom hybridization / lone pairs / angles.
================================================================= */

// 2-letter elements must be matched before 1-letter ones (e.g. "Na" before "N").
const SMILES_ELEMENTS_2 = ["Cl","Br","Li","Na","Mg","Al","Si","Ca","As","Se","Xe",
                            "Fe","Cu","Zn","Ag","Mn","Ni","Co","Cr","Ti","Sn","Pb","Hg","Au","Pt","Pd"];
const SMILES_ELEMENTS_1 = ["B","C","N","O","P","S","F","I","K"];
const KNOWN_ELEMENTS = [...SMILES_ELEMENTS_2, ...SMILES_ELEMENTS_1];

// For transition/post-transition metals, VALENCE_ELECTRONS is set equal to
// TYPICAL_VALENCE (rather than their real d-block electron count) so the
// existing lone-pair formula collapses to 0 lone pairs / steric = bond count —
// a deliberate simplification, not a model of real coordination chemistry
// (VSEPR is already documented as approximate for transition metals).
const VALENCE_ELECTRONS = {
  H:1, B:3, C:4, N:5, O:6, F:7, P:5, S:6, Cl:7, Br:7, I:7,
  Li:1, Na:1, K:1, Mg:2, Ca:2, Al:3, Si:4, As:5, Se:6, Xe:8,
  Fe:3, Cu:2, Zn:2, Ag:1, Mn:2, Ni:2, Co:2, Cr:3, Ti:4, Sn:4, Pb:2, Hg:2, Au:3, Pt:2, Pd:2
};
const TYPICAL_VALENCE = {
  H:1, B:3, C:4, N:3, O:2, F:1, P:3, S:2, Cl:1, Br:1, I:1,
  Li:1, Na:1, K:1, Mg:2, Ca:2, Al:3, Si:4, As:3, Se:2, Xe:0,
  Fe:3, Cu:2, Zn:2, Ag:1, Mn:2, Ni:2, Co:2, Cr:3, Ti:4, Sn:4, Pb:2, Hg:2, Au:3, Pt:2, Pd:2
};
const COLORS = {
  C:"#5cc8ff", O:"#ff7a7a", N:"#7a9bff", S:"#ffd166", P:"#e7a9ff",
  F:"#9fe6ff", Cl:"#9fe6ff", Br:"#9fe6ff", I:"#9fe6ff", H:"#e7ebf2", B:"#ffb070",
  Li:"#ffab91", Na:"#ce93d8", K:"#b39ddb", Mg:"#a5d6a7", Ca:"#81c995", Al:"#cfd8dc",
  Si:"#d7ccc8", As:"#c5cae9", Se:"#ffe082", Xe:"#ea80fc",
  Fe:"#ffb74d", Cu:"#ff8a65", Zn:"#90a4ae", Ag:"#e0e0e0", Mn:"#ab47bc", Ni:"#80cbc4",
  Co:"#64b5f6", Cr:"#4db6ac", Ti:"#b0bec5", Sn:"#bcaaa4", Pb:"#78909c", Hg:"#c6a4d4",
  Au:"#ffd54f", Pt:"#eceff1", Pd:"#b2dfdb"
};

/* ---------------------------------------------------------------
   1. FORMULA PARSER
   Accepts things like "C2H6O", "Ca(OH)2", "NH3", "SF6"
--------------------------------------------------------------- */
function parseFormula(formula){
  formula = formula.trim();
  let i = 0;
  function parseGroup(){
    const counts = {};
    while(i < formula.length && formula[i] !== ')'){
      const c = formula[i];
      if (c === '('){
        i++;
        const inner = parseGroup(); // recurse
        i++; // skip ')'
        let numStr = '';
        while(i<formula.length && /[0-9]/.test(formula[i])){ numStr += formula[i]; i++; }
        const mult = numStr ? parseInt(numStr) : 1;
        for (const el in inner) counts[el] = (counts[el]||0) + inner[el]*mult;
        continue;
      }
      const m = formula.slice(i).match(/^([A-Z][a-z]?)(\d*)/);
      if (!m || !m[1]) { i++; continue; }
      const el = m[1];
      const num = m[2] ? parseInt(m[2]) : 1;
      counts[el] = (counts[el]||0) + num;
      i += m[0].length;
    }
    return counts;
  }
  const counts = parseGroup();
  if (Object.keys(counts).length === 0) throw new Error("No valid elements found.");
  return counts;
}

function canonicalFormula(counts){
  // Hill order: Carbon first, then Hydrogen, then rest alphabetically
  const keys = Object.keys(counts);
  keys.sort((a,b)=>{
    if (a==='C') return -1; if (b==='C') return 1;
    if (a==='H') return -1; if (b==='H') return 1;
    return a.localeCompare(b);
  });
  return keys.map(k => k + (counts[k]>1 ? counts[k] : '')).join('');
}

/* ---------------------------------------------------------------
   2. COMPOUND LIBRARY  (formula -> one or more known structures)
   Each entry stores a SMILES string we already know how to parse.
--------------------------------------------------------------- */
const LIBRARY = [
  { formula:"H2O",      name:"Water",            smiles:"O" },
  { formula:"CH4",      name:"Methane",          smiles:"C" },
  { formula:"NH3",      name:"Ammonia",          smiles:"N" },
  { formula:"CO2",      name:"Carbon dioxide",   smiles:"O=C=O" },
  { formula:"CO",       name:"Carbon monoxide",  smiles:"[C-]#[O+]" },
  { formula:"O3",       name:"Ozone",            smiles:"O=O=O" },
  { formula:"SO2",      name:"Sulfur dioxide",   smiles:"O=S=O" },
  { formula:"SO3",      name:"Sulfur trioxide",  smiles:"O=S(=O)=O" },
  { formula:"H2S",      name:"Hydrogen sulfide", smiles:"S" },
  { formula:"HCl",      name:"Hydrogen chloride",smiles:"Cl" },
  { formula:"HF",       name:"Hydrogen fluoride",smiles:"F" },
  { formula:"BF3",      name:"Boron trifluoride",smiles:"FB(F)F" },
  { formula:"PCl5",     name:"Phosphorus pentachloride", smiles:"ClP(Cl)(Cl)(Cl)Cl" },
  { formula:"SF6",      name:"Sulfur hexafluoride", smiles:"FS(F)(F)(F)(F)F" },
  { formula:"PCl3",     name:"Phosphorus trichloride", smiles:"ClP(Cl)Cl" },
  { formula:"CH2O",     name:"Formaldehyde",     smiles:"C=O" },
  { formula:"C2H2",     name:"Acetylene",        smiles:"C#C" },
  { formula:"C2H4",     name:"Ethylene",         smiles:"C=C" },
  { formula:"C2H6",     name:"Ethane",           smiles:"CC" },
  { formula:"C2H6O",    name:"Ethanol",          smiles:"CCO" },
  { formula:"C2H6O",    name:"Dimethyl ether",   smiles:"COC" },
  { formula:"CH3OH",    name:"Methanol",         smiles:"CO" },
  { formula:"CH4O",     name:"Methanol",         smiles:"CO" },
  { formula:"C2H4O2",   name:"Acetic acid",      smiles:"CC(=O)O" },
  { formula:"C6H6",     name:"Benzene",          smiles:"c1ccccc1" },
  { formula:"C6H12O6",  name:"Glucose",          smiles:"OCC(O)C(O)C(O)C(O)C=O" },
  { formula:"C3H8",     name:"Propane",          smiles:"CCC" },
  { formula:"C3H6O",    name:"Acetone",          smiles:"CC(=O)C" },
  { formula:"C8H18",    name:"Octane",           smiles:"CCCCCCCC" },
  { formula:"N2",       name:"Nitrogen gas",     smiles:"N#N" },
  { formula:"O2",       name:"Oxygen gas",       smiles:"O=O" },
  { formula:"H2O2",     name:"Hydrogen peroxide",smiles:"OO" },
  { formula:"C6H12O6",  name:"Fructose",         smiles:"OCC(=O)C(O)C(O)C(O)CO" },
  { formula:"C3H8O",    name:"Isopropanol",      smiles:"CC(O)C" },
  { formula:"C3H8O",    name:"1-Propanol",       smiles:"CCCO" },
  { formula:"HCN",      name:"Hydrogen cyanide", smiles:"C#N" },
  { formula:"C7H8",     name:"Toluene",          smiles:"Cc1ccccc1" },
];

// Precompute canonical formulas for fast lookup
LIBRARY.forEach(e => { e.canon = canonicalFormula(parseFormula(e.formula)); e.source = 'built-in'; });

// Community-submitted compounds, fetched from the server at load time and
// merged in alongside the built-in LIBRARY for lookup purposes.
let USER_LIBRARY = [];

async function loadUserLibrary(){
  try {
    const res = await fetch('/api/library');
    if (!res.ok) return;
    const data = await res.json();
    USER_LIBRARY = data.map(e => ({ ...e, canon: canonicalFormula(parseFormula(e.formula)), source: e.source === 'pubchem' ? 'pubchem' : 'community' }));
    renderContribTable();
  } catch(e){
    // server not available (e.g. opened as a plain static file) - silently skip
  }
}

function findInLibrary(canon){
  return [...LIBRARY, ...USER_LIBRARY].filter(e => e.canon === canon);
}

/* ---------------------------------------------------------------
   3. SINGLE-CENTRAL-ATOM VSEPR HEURISTIC
   Used when the formula isn't in the library. Works when there's
   exactly ONE non-hydrogen element with a count of 1 (e.g. CO2,
   NH3, SF6, BF3, PCl5) — that becomes the central atom, and every
   other atom attaches to it directly.
--------------------------------------------------------------- */
function tryCentralAtomHeuristic(counts){
  const nonH = Object.keys(counts).filter(e => e !== 'H');
  const candidates = nonH.filter(e => counts[e] === 1);
  if (candidates.length !== 1) return null; // ambiguous or no unique center

  const central = candidates[0];
  // build ligand list: every other non-H atom (each as separate atom) + all H's
  const ligandElements = [];
  nonH.forEach(e=>{
    if (e === central) return;
    for(let k=0;k<counts[e];k++) ligandElements.push(e);
  });
  const hCount = counts['H'] || 0;
  for(let k=0;k<hCount;k++) ligandElements.push('H');

  if (ligandElements.length === 0) return null; // e.g. just "Ne" - nothing to attach

  // Build atoms array directly (skip SMILES, build mol object straight away)
  const atoms = [{ element: central, aromatic:false, charge:0, numH:0 }];
  const bonds = [];
  ligandElements.forEach(el=>{
    atoms.push({ element: el, aromatic:false, charge:0, numH:0, isH: el==='H' });
    bonds.push({ a:0, b:atoms.length-1, order:1, aromatic:false });
  });

  // Distribute extra bond order (double/triple bonds) to satisfy central atom's
  // typical valence, preferring non-hydrogen ligands (H can only form single bonds).
  const centralValence = TYPICAL_VALENCE[central] ?? 4;
  let used = bonds.length; // each ligand currently single-bonded = 1 each
  let deficit = centralValence - used;
  // only spend deficit on non-H ligand bonds, round-robin, max order 3
  let guardCounter = 0;
  while(deficit > 0 && guardCounter < 20){
    guardCounter++;
    let upgraded = false;
    for(const bd of bonds){
      const ligand = atoms[bd.b];
      if (ligand.element === 'H') continue;
      if (bd.order >= 3) continue;
      bd.order += 1;
      deficit -= 1;
      upgraded = true;
      if (deficit <= 0) break;
    }
    if (!upgraded) break; // can't upgrade further (e.g. all H, or maxed out)
  }

  return { atoms, bonds, centralIdeaName: null };
}

/* ---------------------------------------------------------------
   3b. GENERAL CHAIN VSEPR HEURISTIC
   Used when tryCentralAtomHeuristic can't find a single unique center
   (e.g. N2H4, C2H6, S2Cl2 — two or more heavy atoms with no obvious
   single hub). Chains together only the "skeleton-capable" heavy atoms
   (typical valence >= 2, e.g. C/N/O/P/S), since valence-1 atoms like
   halogens can only ever have one neighbour and must sit at the ends as
   ligands, not be chained to each other. Hydrogens and other valence-1
   ligands are then distributed across the chain in valence-capacity
   order (so e.g. propane's middle carbon gets 2 H's, not 3, like its
   end carbons), and chain bond orders are upgraded (single -> double ->
   triple) to satisfy any remaining valence deficit — the same "spend
   the deficit on bond order, not invented atoms" approach as the
   central-atom case.
--------------------------------------------------------------- */
function tryChainHeuristic(counts){
  const nonH = Object.keys(counts).filter(e => e !== 'H');

  const skeletonEls = [];
  const terminalEls = [];
  nonH.forEach(e => {
    const v = TYPICAL_VALENCE[e] ?? 1;
    for (let k=0;k<counts[e];k++){
      (v >= 2 ? skeletonEls : terminalEls).push(e);
    }
  });
  if (skeletonEls.length < 2) return null; // no multi-atom skeleton to chain

  // Higher typical valence first, so more-connected atoms land toward the
  // middle of the chain rather than the ends.
  skeletonEls.sort((a,b) => (TYPICAL_VALENCE[b] ?? 1) - (TYPICAL_VALENCE[a] ?? 1));
  const chainLen = skeletonEls.length;

  const atoms = skeletonEls.map(el => ({ element: el, aromatic:false, charge:0, numH:0 }));
  const bonds = [];
  for (let i=0; i<chainLen-1; i++){
    bonds.push({ a:i, b:i+1, order:1, aromatic:false });
  }

  const usedValence = new Array(chainLen).fill(0);
  bonds.forEach(b => { usedValence[b.a]++; usedValence[b.b]++; });

  // Attach every terminal-only heavy atom (halogens), then every hydrogen,
  // to whichever skeleton atom next has valence room (round-robin, skipping
  // atoms already at capacity).
  const ligands = terminalEls.concat(new Array(counts['H'] || 0).fill('H'));
  let ptr = 0;
  ligands.forEach(el => {
    let idx = -1;
    for (let guard=0; guard<chainLen; guard++){
      const candidate = ptr % chainLen;
      ptr++;
      if (usedValence[candidate] < (TYPICAL_VALENCE[skeletonEls[candidate]] ?? 4)){ idx = candidate; break; }
    }
    if (idx === -1) idx = (ptr - 1) % chainLen; // everyone's full — place anyway rather than drop the atom
    atoms.push({ element: el, aromatic:false, charge:0, numH:0, isH: el==='H' });
    bonds.push({ a: idx, b: atoms.length-1, order:1, aromatic:false });
    usedValence[idx]++;
  });

  // Upgrade chain (skeleton-skeleton) bond orders to satisfy any remaining
  // valence deficit, only where both endpoints still have room.
  let guardCounter = 0, upgraded = true;
  while (upgraded && guardCounter < 30){
    guardCounter++;
    upgraded = false;
    for (let i=0; i<chainLen; i++){
      const cap = TYPICAL_VALENCE[skeletonEls[i]] ?? 4;
      if (usedValence[i] >= cap) continue;
      const bd = bonds.find(b =>
        b.a < chainLen && b.b < chainLen && (b.a === i || b.b === i) && b.order < 3 &&
        usedValence[b.a] < (TYPICAL_VALENCE[skeletonEls[b.a]] ?? 4) &&
        usedValence[b.b] < (TYPICAL_VALENCE[skeletonEls[b.b]] ?? 4)
      );
      if (bd){
        bd.order += 1;
        usedValence[bd.a] += 1;
        usedValence[bd.b] += 1;
        upgraded = true;
      }
    }
  }

  return { atoms, bonds };
}

/* ---------------------------------------------------------------
   4. MINI SMILES PARSER (for library entries)

   Two stages: tokenizeSmiles() turns the raw string into a flat list
   of typed tokens (throwing on anything that isn't valid SMILES
   syntax), then buildMoleculeFromTokens() walks that token list into
   an {atoms, bonds} graph (throwing on invalid *sequencing*, e.g.
   unmatched branches or dangling ring closures). Every place the
   parser used to silently default to a wrong atom or drop a
   character now throws a specific SmilesParseError instead.

   Stereo (@, @@, /, \) and isotopes are parsed so they can't corrupt
   the surrounding element/H-count/charge extraction, then discarded —
   nothing downstream models 3D geometry or isotopic mass, so there's
   nothing to attach that data to.
--------------------------------------------------------------- */
const SMILES_AROMATIC = "bcnops"; // aromatic lowercase organic-subset atoms

class SmilesParseError extends Error {
  constructor(message, pos, smiles){
    super(`${message} (position ${pos} in "${smiles}")`);
    this.name = 'SmilesParseError';
    this.pos = pos;
  }
}

function tokenizeSmiles(smiles){
  const tokens = [];
  let i = 0;
  while (i < smiles.length){
    const c = smiles[i];
    if (c === '('){ tokens.push({ type:'open', pos:i }); i++; continue; }
    if (c === ')'){ tokens.push({ type:'close', pos:i }); i++; continue; }
    if (c === '.'){ tokens.push({ type:'dot', pos:i }); i++; continue; }
    if (c === '-'){ tokens.push({ type:'bond', order:1, aromatic:false, pos:i }); i++; continue; }
    if (c === '='){ tokens.push({ type:'bond', order:2, aromatic:false, pos:i }); i++; continue; }
    if (c === '#'){ tokens.push({ type:'bond', order:3, aromatic:false, pos:i }); i++; continue; }
    if (c === ':'){ tokens.push({ type:'bond', order:1, aromatic:true, pos:i }); i++; continue; }
    if (c === '/' || c === '\\'){ tokens.push({ type:'bondDir', pos:i }); i++; continue; }

    if (c === '%'){
      const digits = smiles.substr(i+1, 2);
      if (!/^\d{2}$/.test(digits)) throw new SmilesParseError(`Expected two digits after '%'`, i, smiles);
      tokens.push({ type:'ring', label: digits, pos:i });
      i += 3; continue;
    }
    if (c >= '0' && c <= '9'){ tokens.push({ type:'ring', label:c, pos:i }); i++; continue; }

    if (c === '['){
      const close = smiles.indexOf(']', i+1);
      if (close === -1) throw new SmilesParseError(`Unclosed '[' bracket atom`, i, smiles);
      tokens.push({ type:'bracketAtom', body: smiles.slice(i+1, close), pos:i });
      i = close + 1; continue;
    }

    if (c === '*'){ tokens.push({ type:'atom', element:'C', aromatic:false, pos:i }); i++; continue; }

    const two = smiles.substr(i, 2);
    if (SMILES_ELEMENTS_2.includes(two)){ tokens.push({ type:'atom', element:two, aromatic:false, pos:i }); i += 2; continue; }
    if (SMILES_ELEMENTS_1.includes(c)){ tokens.push({ type:'atom', element:c, aromatic:false, pos:i }); i++; continue; }
    if (SMILES_AROMATIC.includes(c)){ tokens.push({ type:'atom', element:c.toUpperCase(), aromatic:true, pos:i }); i++; continue; }

    throw new SmilesParseError(`Unrecognized character '${c}'`, i, smiles);
  }
  return tokens;
}

// Bracket-atom grammar (SMILES spec order): isotope? symbol chirality? hcount? charge? (':' class)?
function parseBracketAtom(body, pos, smiles){
  let p = 0;
  const fail = (msg) => { throw new SmilesParseError(`Bad bracket atom '[${body}]': ${msg}`, pos, smiles); };

  // 1. Isotope — leading digits, parsed and discarded (no downstream consumer)
  const isoMatch = body.slice(p).match(/^\d+/);
  if (isoMatch) p += isoMatch[0].length;

  // 2. Symbol
  let element = null, aromatic = false;
  if (body[p] === '*'){ element = 'C'; p += 1; }
  else {
    const two = body.substr(p, 2);
    if (KNOWN_ELEMENTS.includes(two)){ element = two; p += 2; }
    else if (/[A-Z]/.test(body[p] || '')){
      const one = body[p];
      if (!KNOWN_ELEMENTS.includes(one)) fail(`unknown element symbol '${one}'`);
      element = one; p += 1;
    }
    else if (body[p] && SMILES_AROMATIC.includes(body[p])){
      element = body[p].toUpperCase(); aromatic = true; p += 1;
    }
    else fail(`expected an element symbol, got '${body[p] ?? '(end)'}'`);
  }

  // 3. Chirality — parsed and kept on the atom, but nothing downstream reads it
  let chirality = null;
  if (body[p] === '@'){
    if (body[p+1] === '@'){ chirality = '@@'; p += 2; }
    else { chirality = '@'; p += 1; }
  }

  // 4. H-count
  let hCount = 0;
  if (body[p] === 'H'){
    p += 1;
    const digits = body.slice(p).match(/^\d+/);
    if (digits){ hCount = parseInt(digits[0], 10); p += digits[0].length; }
    else hCount = 1;
  }

  // 5. Charge — repeated sign (++, --) or sign+digits (+2) or bare sign
  let charge = 0;
  if (body[p] === '+' || body[p] === '-'){
    const sign = body[p] === '+' ? 1 : -1;
    let q = p + 1;
    if (body[q] === body[p]){
      let count = 1;
      while (body[q] === body[p]){ count++; q++; }
      charge = sign * count; p = q;
    } else {
      const digits = body.slice(q).match(/^\d+/);
      if (digits){ charge = sign * parseInt(digits[0], 10); p = q + digits[0].length; }
      else { charge = sign; p = q; }
    }
  }

  // 6. Atom class — ':' + digits, parsed and discarded
  if (body[p] === ':'){
    p += 1;
    const digits = body.slice(p).match(/^\d+/);
    if (!digits) fail(`expected digits after ':' for atom class`);
    p += digits[0].length;
  }

  if (p !== body.length) fail(`unexpected trailing content '${body.slice(p)}'`);

  return { element, aromatic, charge, chirality, _explicitH: hCount };
}

function buildMoleculeFromTokens(tokens, smiles){
  const atoms = [], bonds = [];
  const ringMap = {}; // label -> { atomIdx, bond }
  const branchStack = [];
  let prev = -1;
  let pendingBond = { order:1, aromatic:false };
  let atomSeen = false;

  // A bond with no explicit symbol (order:1, aromatic:false — the default)
  // between two aromatic atoms is an implicit aromatic bond per SMILES
  // convention (e.g. every ring bond in "c1ccccc1" — nobody writes the
  // explicit ':' symbol in practice), so it should be treated as order 1.5
  // for valence/hydrogen purposes, not as a plain single bond.
  function isDefaultBond(bond){ return bond.order === 1 && !bond.aromatic; }

  function addAtom(atomSpec){
    atoms.push(atomSpec);
    const idx = atoms.length - 1;
    if (prev !== -1){
      const inferAromatic = isDefaultBond(pendingBond) && atoms[prev].aromatic && atomSpec.aromatic;
      bonds.push({ a: prev, b: idx, order: pendingBond.order, aromatic: pendingBond.aromatic || inferAromatic });
    }
    pendingBond = { order:1, aromatic:false };
    prev = idx;
    atomSeen = true;
  }

  function handleRingToken(tok){
    if (prev === -1) throw new SmilesParseError(`Ring bond digit with no preceding atom`, tok.pos, smiles);
    const label = tok.label;
    if (ringMap[label] === undefined){
      ringMap[label] = { atomIdx: prev, bond: pendingBond };
      pendingBond = { order:1, aromatic:false };
    } else {
      const opener = ringMap[label];
      if (opener.atomIdx === prev) throw new SmilesParseError(`Ring bond '${label}' cannot close on the same atom that opened it`, tok.pos, smiles);
      const order = pendingBond.order !== 1 ? pendingBond.order : opener.bond.order;
      const inferAromatic = isDefaultBond(pendingBond) && isDefaultBond(opener.bond) &&
        atoms[opener.atomIdx].aromatic && atoms[prev].aromatic;
      const aromatic = pendingBond.aromatic || opener.bond.aromatic || inferAromatic;
      bonds.push({ a: opener.atomIdx, b: prev, order, aromatic });
      pendingBond = { order:1, aromatic:false };
      delete ringMap[label];
    }
  }

  tokens.forEach(tok => {
    switch (tok.type){
      case 'open':
        if (prev === -1) throw new SmilesParseError(`'(' with no preceding atom to branch from`, tok.pos, smiles);
        branchStack.push(prev);
        break;
      case 'close':
        if (branchStack.length === 0) throw new SmilesParseError(`Unmatched ')' — no open branch to close`, tok.pos, smiles);
        prev = branchStack.pop();
        break;
      case 'dot':
        prev = -1; // start a new disconnected fragment (e.g. ionic salts)
        break;
      case 'bond':
        pendingBond = { order: tok.order, aromatic: tok.aromatic };
        break;
      case 'bondDir':
        break; // parsed so it isn't misread as something else, then discarded
      case 'ring':
        handleRingToken(tok);
        break;
      case 'atom':
        addAtom({ element: tok.element, aromatic: tok.aromatic, charge: 0 });
        break;
      case 'bracketAtom':
        addAtom(parseBracketAtom(tok.body, tok.pos, smiles));
        break;
    }
  });

  if (branchStack.length > 0) throw new SmilesParseError(`Unclosed '(' — missing ')'`, smiles.length, smiles);
  const danglingRing = Object.keys(ringMap)[0];
  if (danglingRing !== undefined) throw new SmilesParseError(`Ring bond label '${danglingRing}' was never closed`, smiles.length, smiles);
  if (!atomSeen) throw new SmilesParseError(`No atoms found in SMILES string`, 0, smiles);

  return { atoms, bonds };
}

function parseSmiles(smiles){
  if (typeof smiles !== 'string' || smiles.trim() === ''){
    throw new SmilesParseError('SMILES string is empty', 0, smiles || '');
  }
  return buildMoleculeFromTokens(tokenizeSmiles(smiles), smiles);
}

/* ---------------------------------------------------------------
   5. ADD IMPLICIT HYDROGENS
--------------------------------------------------------------- */
function addImplicitHydrogens(mol){
  const { atoms, bonds } = mol;
  const bondOrderSum = new Array(atoms.length).fill(0);
  bonds.forEach(b=>{
    const ord = b.aromatic ? 1.5 : b.order;
    bondOrderSum[b.a] += ord;
    bondOrderSum[b.b] += ord;
  });

  for (let idx=0; idx<atoms.length; idx++){
    const atom = atoms[idx];
    if (atom.element === 'H'){ atom.numH = 0; continue; }
    if (atom._explicitH !== undefined){ atom.numH = atom._explicitH; continue; }
    const typical = TYPICAL_VALENCE[atom.element] ?? 4;
    const used = Math.round(bondOrderSum[idx]);
    let h = typical - used - (atom.charge||0);
    atom.numH = Math.max(0, h);
  }

  const newAtoms = atoms.slice();
  const newBonds = bonds.slice();
  atoms.forEach((atom, idx)=>{
    if (atom.element === 'H') return;
    for(let k=0;k<(atom.numH||0);k++){
      newAtoms.push({ element:'H', aromatic:false, charge:0, numH:0, isH:true, parent: idx });
      newBonds.push({ a: idx, b: newAtoms.length-1, order:1, aromatic:false });
    }
  });
  return { atoms:newAtoms, bonds:newBonds };
}

/* ---------------------------------------------------------------
   6. VSEPR ANALYSIS PER HEAVY ATOM
--------------------------------------------------------------- */
function analyze(mol){
  const { atoms, bonds } = mol;
  const neighborBonds = atoms.map(()=>[]);
  bonds.forEach(b=>{ neighborBonds[b.a].push(b); neighborBonds[b.b].push(b); });

  const results = [];
  atoms.forEach((atom, idx)=>{
    if (atom.element === 'H' || atom.isH) return;
    const nb = neighborBonds[idx];
    const sigmaCount = nb.length;
    let totalBondOrder = 0;
    nb.forEach(b=>{ totalBondOrder += (b.aromatic ? 1.5 : b.order); });

    const V = VALENCE_ELECTRONS[atom.element] ?? 4;
    let lonePairs = (V - totalBondOrder - (atom.charge||0)) / 2;
    lonePairs = Math.max(0, Math.round(lonePairs));

    let steric = sigmaCount + lonePairs;
    if (atom.aromatic) steric = Math.max(steric, 3);

    let hybrid, geometry, angles;
    if (atom.aromatic){
      hybrid='sp2'; geometry='trigonal planar (aromatic)'; angles='~120°';
    } else {
      switch(steric){
        case 1: hybrid='s';   geometry='linear (terminal atom)'; angles='—'; break;
        case 2: hybrid='sp';  geometry='linear'; angles='180°'; break;
        case 3: hybrid='sp2'; geometry= lonePairs>0 ? 'bent / trigonal planar e⁻ domain':'trigonal planar'; angles='~120°'; break;
        case 4: hybrid='sp3'; geometry = lonePairs===0?'tetrahedral':lonePairs===1?'trigonal pyramidal':'bent'; angles='~109.5°'; break;
        case 5: hybrid='sp3d'; geometry = lonePairs===0?'trigonal bipyramidal':lonePairs===1?'seesaw':lonePairs===2?'T-shaped':'linear'; angles='90° / 120° / 180°'; break;
        case 6: hybrid='sp3d2'; geometry = lonePairs===0?'octahedral':'square pyramidal / square planar'; angles='90° / 180°'; break;
        default: hybrid='sp3d3+'; geometry='complex/hypervalent'; angles='varies'; break;
      }
    }
    // record neighbor info for the explanation generator
    const neighbors = nb.map(b=>{
      const otherIdx = (b.a === idx) ? b.b : b.a;
      return { element: atoms[otherIdx].element, order: b.aromatic ? '1.5 (aromatic)' : b.order, rawOrder: b.aromatic?1.5:b.order };
    });

    results.push({ idx, element: atom.element, sigmaCount, lonePairs, steric, hybrid, geometry, angles,
                    totalBondOrder, valenceElectrons: V, charge: atom.charge||0, aromatic: atom.aromatic, neighbors });
  });
  return results;
}

/* ---------------------------------------------------------------
   6b. EXPLANATION GENERATOR
   Builds a plain-language, numbers-backed walkthrough of how each
   atom's hybridization/geometry was derived, specific to THIS
   compound's actual bond counts and electron totals.
--------------------------------------------------------------- */
function summarizeNeighbors(neighbors){
  const tally = {};
  neighbors.forEach(n=>{
    const kind = n.order===1?'single':n.order===2?'double':n.order===3?'triple':'aromatic';
    const k = kind+'|'+n.element;
    tally[k] = (tally[k]||0)+1;
  });
  const phrases = [];
  Object.entries(tally).forEach(([k,count])=>{
    const [kind, el] = k.split('|');
    const bondWord = count===1 ? `a ${kind} bond` : `${count} ${kind} bonds`;
    phrases.push(`${bondWord} to ${el}`);
  });
  if (phrases.length === 0) return 'no attached atoms';
  if (phrases.length === 1) return phrases[0];
  return phrases.slice(0,-1).join(', ') + ' and ' + phrases[phrases.length-1];
}

function geometryReason(steric, lonePairs){
  switch(steric){
    case 2: return `A steric number of 2 means the electron domains arrange themselves as far apart as possible in a straight line, giving a linear shape with bond angles of 180°.`;
    case 3: return lonePairs===0
      ? `A steric number of 3 with no lone pairs spreads three bonding domains evenly in a flat triangle (trigonal planar), each ~120° apart.`
      : `A steric number of 3 includes a lone pair taking up one of the three positions in a trigonal-planar electron arrangement, which pushes the remaining bond(s) and gives a bent molecular shape, though the underlying electron-domain angle is still close to 120°.`;
    case 4: return lonePairs===0
      ? `A steric number of 4 with no lone pairs spreads four bonding domains into a tetrahedron, the geometry that minimizes electron repulsion, giving ~109.5° angles.`
      : lonePairs===1
      ? `A steric number of 4 with one lone pair still uses a tetrahedral electron arrangement, but since lone pairs aren't "visible" in the molecular shape, the visible atoms form a trigonal pyramid, with bond angles slightly compressed below the ideal 109.5° by lone-pair repulsion.`
      : `A steric number of 4 with two lone pairs leaves only two bonding positions visible out of a tetrahedral arrangement, producing a bent shape with an angle typically a bit less than 109.5° due to extra lone-pair repulsion.`;
    case 5: return `A steric number of 5 corresponds to a trigonal bipyramidal electron arrangement, with positions at 90°, 120°, and 180° from each other depending on whether they're equatorial or axial; lone pairs preferentially occupy equatorial positions, which is why the resulting shape changes (seesaw, T-shaped, linear) as lone pairs increase.`;
    case 6: return `A steric number of 6 corresponds to an octahedral electron arrangement, with all positions 90° (adjacent) or 180° (opposite) apart; lone pairs occupy positions opposite each other when there are two, giving square-planar geometry, or a single position for square pyramidal.`;
    default: return `This steric number falls outside the simple VSEPR table (1–6 electron domains) and likely reflects an unusual or hypervalent bonding situation.`;
  }
}

function generateExplanations(results){
  return results.map(r=>{
    const neighborText = summarizeNeighbors(r.neighbors);
    const lonePairCalc = `Lone pairs = (valence electrons − total bond order) / 2 = (${r.valenceElectrons} − ${r.totalBondOrder}) / 2 ≈ ${r.lonePairs}.`;
    const stericCalc = `Steric number = σ-bonds (${r.sigmaCount}) + lone pairs (${r.lonePairs}) = ${r.steric}.`;
    const aromaticNote = r.aromatic ? ` This atom is part of an aromatic ring, so it's treated as sp2 regardless of the raw steric-number count.` : '';
    const reason = geometryReason(r.steric, r.lonePairs);
    return `
      <div class="explainBlock">
        <div class="explainHead"><span class="badge ${r.hybrid==='sp'?'sp':r.hybrid==='sp2'?'sp2':r.hybrid==='sp3'?'sp3':r.hybrid==='sp3d'?'sp3d':'sp3d2'}">${r.hybrid}</span> ${r.element}#${r.idx+1}</div>
        <div class="explainBody">
          <p><strong>${r.element}</strong> here has ${neighborText}.
          ${r.element} contributes ${r.valenceElectrons} valence electron(s)${r.charge?` (adjusted for a charge of ${r.charge>0?'+':''}${r.charge})`:''},
          and the total bond order around it (counting each double bond as 2, triple as 3) is ${r.totalBondOrder}.</p>
          <p>${lonePairCalc} ${stericCalc}${aromaticNote}</p>
          <p>${reason}</p>
          <p>That's why this atom is shown as <strong>${r.hybrid}</strong>, with an idealized geometry of <em>${r.geometry}</em> and bond angle(s) of <strong>${r.angles}</strong>.</p>
        </div>
      </div>`;
  }).join('');
}

/* ---------------------------------------------------------------
   7. FORCE-DIRECTED 2D LAYOUT
--------------------------------------------------------------- */
function layout(mol, width=600, height=380){
  const { atoms, bonds } = mol;
  const n = atoms.length;
  const pos = [];
  for(let i=0;i<n;i++){
    pos.push({ x: width/2 + (Math.random()-0.5)*100, y: height/2 + (Math.random()-0.5)*100 });
  }
  const iterations = 400;
  const k = 42;
  for(let it=0; it<iterations; it++){
    const disp = pos.map(()=>({x:0,y:0}));
    for(let a=0;a<n;a++){
      for(let b=a+1;b<n;b++){
        let dx = pos[a].x - pos[b].x, dy = pos[a].y - pos[b].y;
        let dist = Math.sqrt(dx*dx+dy*dy) || 0.01;
        const force = (k*k)/dist;
        const fx = (dx/dist)*force, fy=(dy/dist)*force;
        disp[a].x += fx; disp[a].y += fy;
        disp[b].x -= fx; disp[b].y -= fy;
      }
    }
    bonds.forEach(bd=>{
      let dx = pos[bd.a].x - pos[bd.b].x, dy = pos[bd.a].y - pos[bd.b].y;
      let dist = Math.sqrt(dx*dx+dy*dy) || 0.01;
      const force = (dist*dist)/k;
      const fx = (dx/dist)*force, fy=(dy/dist)*force;
      disp[bd.a].x -= fx; disp[bd.a].y -= fy;
      disp[bd.b].x += fx; disp[bd.b].y += fy;
    });
    const temp = Math.max(1, 12 * (1 - it/iterations));
    for(let i=0;i<n;i++){
      let dx = disp[i].x, dy = disp[i].y;
      let dist = Math.sqrt(dx*dx+dy*dy) || 0.01;
      const capped = Math.min(dist, temp);
      pos[i].x += (dx/dist)*capped;
      pos[i].y += (dy/dist)*capped;
      pos[i].x = Math.min(width-25, Math.max(25, pos[i].x));
      pos[i].y = Math.min(height-25, Math.max(25, pos[i].y));
    }
  }
  return pos;
}

/* ---------------------------------------------------------------
   8. RENDER
--------------------------------------------------------------- */
/* Parse out the first numeric angle value (e.g. "~109.5°" → 109.5, "90° / 120°" → 90) */
function primaryAngle(angleStr){
  const m = angleStr.match(/([\d.]+)°/);
  return m ? parseFloat(m[1]) : null;
}

function render(mol, pos, results){
  const svg = document.getElementById('canvas');
  svg.innerHTML = '';
  const ns = "http://www.w3.org/2000/svg";

  // Build neighbour index: atomIdx -> [atomIdx, ...]
  const neighbours = mol.atoms.map(()=>[]);
  mol.bonds.forEach(b=>{
    neighbours[b.a].push(b.b);
    neighbours[b.b].push(b.a);
  });

  // ---- 1. BOND LINES ----
  mol.bonds.forEach(b=>{
    const p1 = pos[b.a], p2 = pos[b.b];
    const order = b.aromatic ? 1 : b.order;
    const dx = p2.x-p1.x, dy = p2.y-p1.y;
    const len = Math.sqrt(dx*dx+dy*dy) || 1;
    const nx = -dy/len, ny = dx/len;
    const offsets = order===1 ? [0] : order===2 ? [-3,3] : [-5,0,5];
    offsets.forEach(off=>{
      const line = document.createElementNS(ns,"line");
      line.setAttribute("x1", p1.x + nx*off);
      line.setAttribute("y1", p1.y + ny*off);
      line.setAttribute("x2", p2.x + nx*off);
      line.setAttribute("y2", p2.y + ny*off);
      line.setAttribute("stroke", b.aromatic ? "#888" : "#aab2c5");
      line.setAttribute("stroke-width", "2");
      line.setAttribute("stroke-dasharray", b.aromatic ? "4,3" : "none");
      svg.appendChild(line);
    });
  });

  // ---- 2. BOND-ANGLE ARCS ----
  // Build a lookup from atom index -> VSEPR result (only for non-H atoms with ≥2 neighbours)
  const resultByIdx = {};
  if (results){
    results.forEach(r => { resultByIdx[r.idx] = r; });
  }

  mol.atoms.forEach((atom, idx)=>{
    if (atom.element === 'H' || atom.isH) return;
    const nb = neighbours[idx];
    // Need at least 2 neighbours for an angle; ignore lone-atom terminal (H, halogens at end)
    if (nb.length < 2) return;
    const res = resultByIdx[idx];
    if (!res) return;
    const idealAngleStr = res.angles;
    if (idealAngleStr === '—') return;

    const cx = pos[idx].x, cy = pos[idx].y;

    // Sort neighbours by their 2D angle around this atom so arcs are drawn
    // between adjacent bond pairs (cleanest visual result).
    const sorted = nb
      .map(ni => ({ ni, a: Math.atan2(pos[ni].y - cy, pos[ni].x - cx) }))
      .sort((x,y) => x.a - y.a);

    // For each adjacent pair in the sorted ring, draw one arc + label.
    // When multiple bond pairs share the same ideal angle value (e.g. every
    // ~109.5° pair in a tetrahedral center), only the first is drawn — the
    // rest are visually redundant.
    const drawnAngleLabels = new Set();

    sorted.forEach((cur, i)=>{
      const next = sorted[(i+1) % sorted.length];
      const a1 = cur.a, a2 = next.a;

      // Angular span going counter-clockwise from a1 to a2
      let span = a2 - a1;
      if (span <= 0) span += 2*Math.PI;

      // Skip very wide gaps (>270°): they're the "back" of the arc and
      // would produce misleading large arcs on a 2D projection.
      if (span > 4.71) return; // 270° in radians

      // Determine the ideal angle label for THIS specific arc span.
      // For atoms with a single ideal angle value, use that for all arcs.
      // For multi-angle atoms (sp3d, sp3d2), pick the closest ideal angle to the drawn span.
      let angleLabel = idealAngleStr;
      const spanDeg = span * 180 / Math.PI;
      const allAngles = idealAngleStr.match(/([\d.]+)°/g);
      if (allAngles && allAngles.length > 1){
        let best = allAngles[0];
        let bestDiff = Infinity;
        allAngles.forEach(a=>{
          const val = parseFloat(a);
          const diff = Math.abs(val - spanDeg);
          if (diff < bestDiff){ bestDiff = diff; best = a; }
        });
        angleLabel = best;
      } else if (allAngles && allAngles.length === 1){
        angleLabel = allAngles[0];
      }

      // Skip this arc entirely (line + label) if this ideal angle value has
      // already been drawn once for this atom.
      if (drawnAngleLabels.has(angleLabel)) return;
      drawnAngleLabels.add(angleLabel);

      // Pick a radius that sits just outside the atom circle
      const arcR = 22;
      const sx = cx + arcR * Math.cos(a1);
      const sy = cy + arcR * Math.sin(a1);
      const ex = cx + arcR * Math.cos(a2);
      const ey = cy + arcR * Math.sin(a2);
      const largeArc = span > Math.PI ? 1 : 0;
      // sweep=1 means clockwise / increasing angle direction
      const sweep = 1;

      const arc = document.createElementNS(ns, "path");
      arc.setAttribute("d", `M ${sx} ${sy} A ${arcR} ${arcR} 0 ${largeArc} ${sweep} ${ex} ${ey}`);
      arc.setAttribute("fill","none");
      arc.setAttribute("stroke","#ffd166");
      arc.setAttribute("stroke-width","1.2");
      arc.setAttribute("stroke-dasharray","3,2");
      svg.appendChild(arc);

      // Label: place at midpoint angle, slightly further out
      const midAngle = a1 + span/2;
      const labelR = arcR + 13;
      const lx = cx + labelR * Math.cos(midAngle);
      const ly = cy + labelR * Math.sin(midAngle);

      const text = document.createElementNS(ns, "text");
      text.setAttribute("x", lx);
      text.setAttribute("y", ly + 3.5);
      text.setAttribute("text-anchor","middle");
      text.setAttribute("font-size","9");
      text.setAttribute("font-weight","700");
      text.setAttribute("fill","#ffd166");
      // Thin white halo for legibility over dark bonds
      text.setAttribute("paint-order","stroke");
      text.setAttribute("stroke","#0f1115");
      text.setAttribute("stroke-width","2.5");
      text.setAttribute("stroke-linejoin","round");
      text.textContent = angleLabel;
      svg.appendChild(text);
    });
  });

  // ---- 3. ATOM CIRCLES + LABELS ----
  mol.atoms.forEach((atom, idx)=>{
    const p = pos[idx];
    const g = document.createElementNS(ns,"g");
    const r = (atom.element==='H') ? 8 : 14;
    const circle = document.createElementNS(ns,"circle");
    circle.setAttribute("cx", p.x); circle.setAttribute("cy", p.y); circle.setAttribute("r", r);
    circle.setAttribute("fill", COLORS[atom.element] || "#888");
    circle.setAttribute("stroke", "#0f1115"); circle.setAttribute("stroke-width","1.5");
    g.appendChild(circle);
    const label = document.createElementNS(ns,"text");
    label.setAttribute("x", p.x); label.setAttribute("y", p.y+4);
    label.setAttribute("text-anchor","middle");
    label.setAttribute("font-size", atom.element==='H' ? "9" : "11");
    label.setAttribute("font-weight","700");
    label.setAttribute("fill", "#0b0d12");
    label.textContent = atom.element;
    g.appendChild(label);
    svg.appendChild(g);
  });
}

/* ---------------------------------------------------------------
   9. DRIVER
--------------------------------------------------------------- */
function badgeClass(h){
  if (h==='sp') return 'sp';
  if (h==='sp2') return 'sp2';
  if (h==='sp3') return 'sp3';
  if (h==='sp3d') return 'sp3d';
  if (h==='sp3d2') return 'sp3d2';
  return 'sp3';
}

let lastMol = null, lastResults = null;

function buildAndShow(mol, label){
  mol = addImplicitHydrogens(mol);
  const pos = layout(mol);
  const results = analyze(mol);
  render(mol, pos, results);
  lastMol = mol;
  lastResults = results;
  const tbody = document.querySelector('#resultTable tbody');
  tbody.innerHTML = '';
  results.forEach(r=>{
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${r.element}<span style="color:var(--muted)">#${r.idx+1}</span></td>
      <td>${r.sigmaCount}</td>
      <td>${r.lonePairs}</td>
      <td>${r.steric}</td>
      <td><span class="badge ${badgeClass(r.hybrid)}">${r.hybrid}</span></td>
      <td>${r.geometry} · <strong>${r.angles}</strong></td>
    `;
    tbody.appendChild(tr);
  });
  // A null label means the caller (e.g. the isomer picker) already populated
  // isomerNote itself and is just asking us to render the default structure —
  // don't wipe out what it just built.
  if (label) document.getElementById('isomerNote').textContent = label;

  document.getElementById('explainWrap').innerHTML = generateExplanations(results);
}

function showSingleMatch(match, note){
  const mol = parseSmiles(match.smiles);
  buildAndShow(mol, note);
}

function showIsomerPicker(matches, canon, introHtml){
  const note = document.getElementById('isomerNote');
  note.innerHTML = introHtml;
  matches.forEach(m=>{
    const btn = document.createElement('button');
    btn.className = 'secondary isomerBtn';
    btn.textContent = m.name;
    btn.onclick = ()=>{
      const mol = parseSmiles(m.smiles);
      buildAndShow(mol, `Showing ${m.name} (${canon}).`);
    };
    note.appendChild(btn);
  });
  // show the first one by default
  const mol = parseSmiles(matches[0].smiles);
  buildAndShow(mol, null);
}

// Ask the server (which proxies PubChem) for compounds matching this
// canonical formula. Returns [] on any failure (offline, PubChem down, no hit).
async function lookupOnline(canon){
  try {
    const res = await fetch('/api/lookup/' + encodeURIComponent(canon));
    if (!res.ok) return [];
    const data = await res.json();
    const compounds = (data.compounds || []).map(c => ({ name: c.name, smiles: c.smiles }));
    // PubChem's formula search also returns isotopologues (deuterated,
    // tritiated, etc.) that share the same connectivity SMILES as the
    // parent compound — collapse those down to one entry.
    const seenSmiles = new Set();
    return compounds.filter(c => {
      if (seenSmiles.has(c.smiles)) return false;
      seenSmiles.add(c.smiles);
      return true;
    });
  } catch(e){
    return [];
  }
}

// Persist newly-discovered compounds into the shared library so future
// lookups for this formula resolve locally instead of hitting PubChem again.
async function cacheWebMatches(matches, canon){
  for (const m of matches){
    try {
      await fetch('/api/library', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name: m.name, formula: canon, smiles: m.smiles, source: 'pubchem' })
      });
    } catch(e){ /* best-effort cache; a failed save just means we ask PubChem again next time */ }
  }
  await loadUserLibrary();
}

async function run(){
  const errEl = document.getElementById('err');
  errEl.textContent = '';
  document.getElementById('isomerNote').textContent = '';
  const formula = document.getElementById('formula').value.trim();
  if (!formula) return;

  try {
    const counts = parseFormula(formula);
    const canon = canonicalFormula(counts);
    const matches = findInLibrary(canon);

    if (matches.length === 1){
      showSingleMatch(matches[0], `Recognized as ${matches[0].name} (${canon}).`);
      return;
    }

    if (matches.length > 1){
      showIsomerPicker(matches, canon, `Multiple known compounds share the formula <strong>${canon}</strong> (structural isomers). Pick one: `);
      return;
    }

    // Not in the local library -> check PubChem for the real structure
    // before falling back to our own VSEPR/valence-rule guesses.
    errEl.textContent = `"${canon}" isn't in the local library — checking PubChem...`;
    const webMatches = await lookupOnline(canon);
    if (webMatches.length > 0){
      errEl.textContent = '';
      await cacheWebMatches(webMatches, canon);
      if (webMatches.length === 1){
        showSingleMatch(webMatches[0], `Not in the local library — found "${webMatches[0].name}" (${canon}) via PubChem and cached for next time.`);
      } else {
        showIsomerPicker(webMatches, canon, `Not in the local library — PubChem found multiple compounds matching <strong>${canon}</strong> (now cached for next time). Pick one: `);
      }
      return;
    }
    errEl.textContent = '';

    // Not on PubChem either -> fall back to our own VSEPR/valence-rule
    // structural estimate. Try a single unique central atom first (precise
    // when it applies, e.g. CO2, NH3, SF6), but don't stop there — if there's
    // no unique center, build a general multi-atom chain instead (e.g. N2H4,
    // C2H6, S2Cl2) rather than giving up.
    const heuristicMol = tryCentralAtomHeuristic(counts);
    if (heuristicMol){
      buildAndShow(heuristicMol, `No match in the library or on PubChem — structure estimated from VSEPR/valence rules (single-center) for ${canon}.`);
      return;
    }

    const chainMol = tryChainHeuristic(counts);
    if (chainMol){
      buildAndShow(chainMol, `No match in the library or on PubChem — structure estimated from VSEPR/valence rules (chain) for ${canon}.`);
      return;
    }

    throw new Error(
      `"${canon}" isn't in the built-in compound library, wasn't found on PubChem, and doesn't fit ` +
      `a central-atom or chain structure that VSEPR/valence rules can resolve. Try a SMILES-capable ` +
      `lookup, or pick a known compound from the presets.`
    );

  } catch(e){
    errEl.textContent = "Couldn't analyze that formula: " + e.message;
    document.querySelector('#resultTable tbody').innerHTML = '';
    document.getElementById('canvas').innerHTML = '';
    document.getElementById('explainWrap').innerHTML = '';
  }
}

/* ---------------------------------------------------------------
   9b. 3D STRUCTURE VIEW
   Generates 3D coordinates from the same idealized VSEPR geometry
   already computed per atom (steric number -> a fixed set of
   idealized bond directions), not a physically minimized structure —
   consistent with the 2D diagram's own "idealized angle" labels.
   Placement is a simple BFS: each atom's direction template is
   rotated so one slot points back at its parent, and the remaining
   slots go to its not-yet-placed neighbours. Ring-closing bonds (to
   an already-placed atom) are just drawn, not used to move anything.
   Rendered via 3Dmol.js as a MOL (V2000) block built from the same
   atoms/bonds this app already models.
--------------------------------------------------------------- */
const VSEPR_DIRECTIONS = {
  1: [[0,0,1]],
  2: [[0,0,1],[0,0,-1]],
  3: [[1,0,0],[-0.5,0.8660254,0],[-0.5,-0.8660254,0]],
  4: [[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],
  5: [[0,0,1],[0,0,-1],[1,0,0],[-0.5,0.8660254,0],[-0.5,-0.8660254,0]],
  6: [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],
};

function vNormalize(v){
  const len = Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) || 1;
  return [v[0]/len, v[1]/len, v[2]/len];
}
function vAdd(a,b){ return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
function vScale(a,s){ return [a[0]*s, a[1]*s, a[2]*s]; }
function vCross(a,b){ return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]; }
function vDot(a,b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

// Returns a function that rotates any vector the same way `from` (unit) needs
// to rotate to land on `to` (unit) — Rodrigues' rotation formula.
function rotationAligning(from, to){
  const f = vNormalize(from), t = vNormalize(to);
  const c = vDot(f, t);
  if (c > 0.9999) return (v) => v;
  if (c < -0.9999){
    let axis = vCross(f, [1,0,0]);
    if (vDot(axis,axis) < 1e-6) axis = vCross(f, [0,1,0]);
    axis = vNormalize(axis);
    return (v) => vAdd(vScale(axis, 2*vDot(axis, v)), vScale(v, -1));
  }
  const axis = vNormalize(vCross(f, t));
  const s = Math.sqrt(1 - c*c);
  return (v) => vAdd(vAdd(vScale(v, c), vScale(vCross(axis, v), s)), vScale(axis, vDot(axis, v) * (1 - c)));
}

function layout3D(mol){
  const { atoms, bonds } = mol;
  const neighbours = atoms.map(() => []);
  bonds.forEach(b => { neighbours[b.a].push(b.b); neighbours[b.b].push(b.a); });

  const resultByIdx = {};
  (lastResults || []).forEach(r => { resultByIdx[r.idx] = r; });

  const pos = new Array(atoms.length).fill(null);
  const visited = new Array(atoms.length).fill(false);
  const BOND_LEN = 1.5;

  function stericOf(idx){
    const r = resultByIdx[idx];
    return r ? Math.min(Math.max(r.steric, 1), 6) : 1; // terminal atoms (H, etc.) just need one slot
  }

  for (let start = 0; start < atoms.length; start++){
    if (visited[start]) continue;
    pos[start] = [0,0,0];
    visited[start] = true;
    const queue = [{ idx: start, parent: -1, inDir: [0,0,1] }];
    while (queue.length){
      const { idx, parent, inDir } = queue.shift();
      const template = VSEPR_DIRECTIONS[stericOf(idx)] || VSEPR_DIRECTIONS[4];
      const backToParent = parent === -1 ? [0,0,1] : vNormalize(vScale(inDir, -1));
      const rotate = rotationAligning(template[0], backToParent);
      const rotatedTemplate = template.map(rotate);

      // Slot 0 is reserved for the bond back to the parent — except for the
      // root atom, which has no parent bond and so gets every slot free.
      let slot = parent === -1 ? 0 : 1;
      neighbours[idx].forEach(n => {
        if (n === parent || visited[n]) return; // already placed (ring closure) - just a bond, no new position
        const dir = rotatedTemplate[slot] || rotatedTemplate[rotatedTemplate.length-1] || [0,0,1];
        slot++;
        pos[n] = vAdd(pos[idx], vScale(dir, BOND_LEN));
        visited[n] = true;
        queue.push({ idx: n, parent: idx, inDir: dir });
      });
    }
  }
  return pos;
}

function buildMolBlock(mol, pos){
  const { atoms, bonds } = mol;
  const pad = (s, n) => String(s).padStart(n);
  const lines = ['', '  StructExpl', '', `${pad(atoms.length,3)}${pad(bonds.length,3)}  0  0  0  0  0  0  0  0999 V2000`];
  atoms.forEach((a, i) => {
    const p = pos[i] || [0,0,0];
    lines.push(`${p[0].toFixed(4).padStart(10)}${p[1].toFixed(4).padStart(10)}${p[2].toFixed(4).padStart(10)} ${a.element.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
  });
  bonds.forEach(b => {
    const order = Math.max(1, Math.min(3, Math.round(b.order))); // MOL V2000 has no 1.5 (aromatic) order
    lines.push(`${pad(b.a+1,3)}${pad(b.b+1,3)}${pad(order,3)}  0`);
  });
  lines.push('M  END');
  return lines.join('\n');
}

let viewer3d = null;

function show3D(){
  if (!lastMol){ return; }
  const overlay = document.getElementById('viewer3dOverlay');
  const container = document.getElementById('viewer3dContainer');
  overlay.style.display = 'flex';
  container.innerHTML = '';

  const pos = layout3D(lastMol);
  const molBlock = buildMolBlock(lastMol, pos);
  const elementColor = (atom) => COLORS[atom.elem] || '#888888';

  viewer3d = $3Dmol.createViewer(container, { backgroundColor: '#171a21' });
  viewer3d.addModel(molBlock, 'sdf');
  viewer3d.setStyle({}, {
    stick: { radius: 0.14, colorfunc: elementColor },
    sphere: { scale: 0.28, colorfunc: elementColor }
  });
  viewer3d.zoomTo();
  viewer3d.render();

  render3DKey(lastMol);
}

function render3DKey(mol){
  const keyEl = document.getElementById('viewer3dKey');
  const seen = new Set();
  const elements = [];
  mol.atoms.forEach(a => { if (!seen.has(a.element)){ seen.add(a.element); elements.push(a.element); } });
  keyEl.innerHTML = elements.map(el =>
    `<span><span class="dot" style="background:${COLORS[el] || '#888888'}"></span>${el}</span>`
  ).join('');
}

function close3D(){
  document.getElementById('viewer3dOverlay').style.display = 'none';
  viewer3d = null;
}

/* ---------------------------------------------------------------
   10. PRESET BUTTONS
--------------------------------------------------------------- */
const PRESETS = [
  ["Water","H2O"], ["Methane","CH4"], ["Ammonia","NH3"], ["CO2","CO2"],
  ["SO2","SO2"], ["SF6","SF6"], ["PCl5","PCl5"], ["BF3","BF3"],
  ["Ethanol/DME","C2H6O"], ["Benzene","C6H6"], ["Acetic acid","C2H4O2"],
  ["Glucose","C6H12O6"], ["Ethylene","C2H4"], ["Acetylene","C2H2"],
];
const presetsDiv = document.getElementById('presets');
PRESETS.forEach(([name, formula])=>{
  const btn = document.createElement('button');
  btn.className = 'secondary';
  btn.textContent = name;
  btn.onclick = ()=>{ document.getElementById('formula').value = formula; run(); };
  presetsDiv.appendChild(btn);
});

document.getElementById('formula').addEventListener('keydown', (e)=>{
  if (e.key === 'Enter') run();
});

/* ---------------------------------------------------------------
   11. TAB SWITCHING
--------------------------------------------------------------- */
document.querySelectorAll('.tabBtn').forEach(btn=>{
  btn.addEventListener('click', ()=>{
    document.querySelectorAll('.tabBtn').forEach(b=>b.classList.remove('active'));
    document.querySelectorAll('.tabPanel').forEach(p=>p.style.display='none');
    btn.classList.add('active');
    document.getElementById('tab-' + btn.dataset.tab).style.display = '';
  });
});

/* ---------------------------------------------------------------
   12. CONTRIBUTE-A-COMPOUND FORM
--------------------------------------------------------------- */
function computeFormulaFromSmiles(smiles){
  const mol = addImplicitHydrogens(parseSmiles(smiles));
  const counts = {};
  mol.atoms.forEach(a => { counts[a.element] = (counts[a.element]||0) + 1; });
  return canonicalFormula(counts);
}

function renderContribTable(){
  const tbody = document.querySelector('#contribTable tbody');
  if (!tbody) return;
  const all = [...LIBRARY, ...USER_LIBRARY];
  tbody.innerHTML = all.map(e => `
    <tr>
      <td>${e.name}</td>
      <td>${e.formula}</td>
      <td style="font-family:Consolas,monospace;">${e.smiles}</td>
      <td>${e.source === 'pubchem' ? '🌐 PubChem' : e.source === 'community' ? '🧪 community' : 'built-in'}</td>
    </tr>`).join('');
}

async function submitContribution(){
  const msgEl = document.getElementById('contribMsg');
  msgEl.className = 'note';
  msgEl.textContent = '';

  const name = document.getElementById('newName').value.trim();
  const formulaInput = document.getElementById('newFormula').value.trim();
  const smiles = document.getElementById('newSmiles').value.trim();

  if (!name || !formulaInput || !smiles){
    msgEl.className = 'note error';
    msgEl.textContent = 'Please fill in a name, formula, and SMILES string.';
    return;
  }

  let typedCanon, derivedCanon;
  try {
    typedCanon = canonicalFormula(parseFormula(formulaInput));
  } catch(e){
    msgEl.className = 'note error';
    msgEl.textContent = "Couldn't parse that formula: " + e.message;
    return;
  }
  try {
    derivedCanon = computeFormulaFromSmiles(smiles);
  } catch(e){
    msgEl.className = 'note error';
    msgEl.textContent = "Couldn't parse that SMILES string: " + e.message;
    return;
  }

  if (typedCanon !== derivedCanon){
    msgEl.className = 'note error';
    msgEl.textContent = `That SMILES implies the formula ${derivedCanon}, which doesn't match the formula you typed (${typedCanon}). Double check both fields.`;
    return;
  }

  try {
    const res = await fetch('/api/library', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, formula: formulaInput, smiles })
    });
    if (!res.ok){
      const err = await res.json().catch(()=>({error:'Unknown error'}));
      throw new Error(err.error || 'Server rejected the submission.');
    }
    msgEl.className = 'note successMsg';
    msgEl.textContent = `Added "${name}" (${typedCanon}) to the shared library. Anyone using this app can now type "${formulaInput}" to find it.`;
    document.getElementById('newName').value = '';
    document.getElementById('newFormula').value = '';
    document.getElementById('newSmiles').value = '';
    await loadUserLibrary();
  } catch(e){
    msgEl.className = 'note error';
    msgEl.textContent = "Couldn't save to the server: " + e.message +
      " (Note: this feature needs the Node server from server.js running — it won't work if this page was opened as a plain static file.)";
  }
}

/* ---------------------------------------------------------------
   13. INIT
--------------------------------------------------------------- */
loadUserLibrary();
run();

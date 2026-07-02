const express = require('express');
const fs = require('fs');
const path = require('path');

const app = express();
app.use(express.json());
app.use(express.static(path.join(__dirname, 'public')));

const LIB_FILE = path.join(__dirname, 'library.json');

function readLibrary(){
  try {
    return JSON.parse(fs.readFileSync(LIB_FILE, 'utf8'));
  } catch (e) {
    return [];
  }
}

function writeLibrary(data){
  fs.writeFileSync(LIB_FILE, JSON.stringify(data, null, 2));
}

// GET the full community-submitted library
app.get('/api/library', (req, res) => {
  res.json(readLibrary());
});

// POST a new compound: { name, formula, smiles, source? }
app.post('/api/library', (req, res) => {
  const { name, formula, smiles, source } = req.body || {};
  if (!name || !formula || !smiles){
    return res.status(400).json({ error: 'name, formula, and smiles are all required.' });
  }
  if (String(name).length > 100 || String(formula).length > 60 || String(smiles).length > 300){
    return res.status(400).json({ error: 'One of the fields is too long.' });
  }
  const library = readLibrary();
  library.push({
    name: String(name).trim(),
    formula: String(formula).trim(),
    smiles: String(smiles).trim(),
    source: source === 'pubchem' ? 'pubchem' : 'community',
    addedAt: new Date().toISOString()
  });
  writeLibrary(library);
  res.json({ ok: true, count: library.length });
});

// GET compounds matching a molecular formula from PubChem — used as a
// fallback when a formula isn't in the built-in or community library.
app.get('/api/lookup/:formula', async (req, res) => {
  const formula = String(req.params.formula || '').trim();
  if (!formula || formula.length > 60){
    return res.status(400).json({ error: 'Invalid formula.', compounds: [] });
  }
  try {
    const cidRes = await fetch(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/${encodeURIComponent(formula)}/cids/JSON`
    );
    if (!cidRes.ok) return res.json({ compounds: [] });
    const cidData = await cidRes.json();
    const cids = (cidData.IdentifierList && cidData.IdentifierList.CID || []).slice(0, 5);
    if (cids.length === 0) return res.json({ compounds: [] });

    const propRes = await fetch(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cids.join(',')}/property/IUPACName,ConnectivitySMILES/JSON`
    );
    if (!propRes.ok) return res.json({ compounds: [] });
    const propData = await propRes.json();
    const properties = (propData.PropertyTable && propData.PropertyTable.Properties) || [];
    const compounds = properties
      .filter(p => p.ConnectivitySMILES && p.IUPACName)
      .map(p => ({ name: p.IUPACName, smiles: p.ConnectivitySMILES, cid: p.CID }));
    res.json({ compounds });
  } catch (e) {
    res.status(502).json({ error: 'Could not reach PubChem.', compounds: [] });
  }
});

const PORT = process.env.PORT || 8080;
app.listen(PORT, () => {
  console.log(`Molecule explorer running on port ${PORT}`);
});

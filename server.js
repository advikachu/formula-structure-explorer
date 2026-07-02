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

// POST a new compound: { name, formula, smiles }
app.post('/api/library', (req, res) => {
  const { name, formula, smiles } = req.body || {};
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
    addedAt: new Date().toISOString()
  });
  writeLibrary(library);
  res.json({ ok: true, count: library.length });
});

const PORT = process.env.PORT || 8080;
app.listen(PORT, () => {
  console.log(`Molecule explorer running on port ${PORT}`);
});

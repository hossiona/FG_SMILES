#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import io
import base64
import pandas as pd
from flask import Flask, render_template, request, Response
from rdkit import Chem
from rdkit.Chem import Draw

app = Flask(__name__)

# --- DEFINE PATTERNS & FUNCTIONS ---
functional_groups = {
    "Carboxylic Acid": "[CX3](=O)[OX2H1]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Aldehyde": "[CX3H](=O)[#6]",
    "Alcohol": "[#6][OX2H1]",
    "Phenol": "[OX2H][c1ccccc1]",
    "Ether": "[#6][OX2][#6]",
    "Sulfide": "[#6][SX2][#6]",
    "Alkene": "C=C",
    "Alkyne": "C#C",
    "Ester": "[CX3](=O)[OX2H0][#6]",
    "Amine": "[NX3;H2,H1,H0;!$(NC=O)]",
    "Amide": "[CX3](=O)[NX3]",
    "Benzene": "c1ccccc1",
    "Anhydride": "[CX3](=O)[OX2][CX3](=O)",
    "AlkylHalide": "[#6][F,Cl,Br,I]",
    "Cycloalkene": "[C;R]=[C;R]",
    "Thiol": "[#6][SX2H]",
    "Disulfide": "[#6][SX2][SX2][#6]"
}

def identify_functional_groups(mol):
    if mol is None:
        return "Invalid SMILES"

    found_groups = []
    for group_name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            found_groups.append(f"{group_name} ({len(matches)})")

    return ", ".join(found_groups) if found_groups else "None Detected"

def mol_to_base64(mol):
    """Converts an RDKit molecule into a base64 encoded PNG string."""
    img = Draw.MolToImage(mol, size=(400, 400))
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
    return f"data:image/png;base64,{img_str}"

# --- ROUTES ---

@app.route("/", methods=["GET", "POST"])
def index():
    # Setup default values
    smiles_input = "CC(=O)Oc1ccccc1C(=O)O"
    image_data = None
    groups_list = None
    error_msg = None

    if request.method == "POST" and "analyze_single" in request.form:
        smiles_input = request.form.get("smiles_input", "").strip()
        mol = Chem.MolFromSmiles(smiles_input)

        if mol is None:
            error_msg = "❌ Invalid SMILES string. RDKit could not parse the structure."
        else:
            image_data = mol_to_base64(mol)
            groups_str = identify_functional_groups(mol)
            groups_list = groups_str.split(", ") if groups_str != "None Detected" else []

    return render_template("index.html", smiles_input=smiles_input, image_data=image_data, groups_list=groups_list, error_msg=error_msg)

@app.route("/batch", methods=["POST"])
def batch_processing():
    if "csv_file" not in request.files:
        return "No file part", 400

    file = request.files["csv_file"]
    if file.filename == "":
        return "No selected file", 400

    if file and file.filename.endswith(".csv"):
        df = pd.read_csv(file)
        smiles_col = request.form.get("smiles_col")

        if smiles_col not in df.columns:
            return f"Column '{smiles_col}' not found in CSV.", 400

        # Run chemical batch calculation
        mols = df[smiles_col].apply(Chem.MolFromSmiles)
        df['Detected Functional Groups'] = mols.apply(identify_functional_groups)

        # Stream the CSV directly back down as an attachment download
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)

        return Response(
            csv_buffer.getvalue(),
            mimetype="text/csv",
            headers={"Content-disposition": "attachment; filename=analyzed_molecules.csv"}
        )

    return "Invalid file format", 400

if __name__ == "__main__":
    app.run(debug=True)


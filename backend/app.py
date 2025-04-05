from flask import Flask, request, jsonify, render_template
from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib
import pandas as pd

app = Flask(__name__, template_folder="templates", static_folder="static")

# Load trained model (for example purposes, replace with your actual model)
model = joblib.load("model.joblib")

# Mock dictionary for molecule generation
MOCK_MOLECULES = {
    "cancer": "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin
    "diabetes": "CN1C=NC2=C1C(=O)N(C(=O)N2)C",  # Metformin
    "fever": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",   # Ibuprofen
}

def extract_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {
        "Molecular_Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "H_Bond_Donors": Descriptors.NumHDonors(mol),
        "H_Bond_Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol)
    }

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/generate_molecule", methods=["POST"])
def generate_molecule():
    data = request.get_json()
    disease = data.get("disease_target", "").lower()

    smiles = MOCK_MOLECULES.get(disease)
    if not smiles:
        return jsonify({"error": "No molecule found for this disease."}), 400

    return jsonify({"molecule": smiles})

@app.route("/analyze_molecule", methods=["POST"])
def analyze_molecule():
    data = request.get_json()
    smiles = data.get("molecule")
    features = extract_features(smiles)

    if features is None:
        return jsonify({"error": "Invalid SMILES input."}), 400

    df = pd.DataFrame([features])
    prediction = model.predict(df)[0]  # Binary prediction (e.g., 1 for likely drug)
    return jsonify({"likelihood": "Likely Drug" if prediction else "Unlikely Drug"})

@app.route("/optimize_molecule", methods=["POST"])
def optimize_molecule():
    data = request.get_json()
    smiles = data.get("molecule")

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return jsonify({"error": "Invalid molecule."}), 400

    # Simple "optimization": remove H atoms to reduce MW
    try:
        mol_noh = Chem.RemoveHs(mol)
        new_smiles = Chem.MolToSmiles(mol_noh)
        return jsonify({"optimized": new_smiles})
    except Exception as e:
        return jsonify({"error": f"Optimization failed: {str(e)}"}), 500

@app.route("/generate_report", methods=["POST"])
def generate_report():
    data = request.get_json()
    smiles = data.get("molecule")
    features = extract_features(smiles)

    if features is None:
        return jsonify({"error": "Invalid molecule."}), 400

    # Generate a simple report as string
    report = "\n".join([f"{k}: {v:.2f}" for k, v in features.items()])
    return jsonify({"report": report})

if __name__ == "__main__":
    app.run(debug=True)

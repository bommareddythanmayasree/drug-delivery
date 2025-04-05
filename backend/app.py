from flask import Flask, render_template, request, jsonify
from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib
import pandas as pd

app = Flask(__name__, template_folder="templates", static_folder="static")

# Load the trained ML model
model = joblib.load("model.joblib")

def get_molecular_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {
        "Molecular_Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "H_Bond_Donors": Descriptors.NumHDonors(mol),
        "H_Bond_Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
    }

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/generate_molecule", methods=["POST"])
def generate_molecule():
    data = request.get_json()
    disease_name = data.get("disease", "").lower()

    # Example SMILES mapping (mock logic)
    mock_molecules = {
        "cancer": "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin
        "diabetes": "CN1C=NC2=C1C(=O)N(C(=O)N2)C",  # Metformin
        "fever": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",   # Ibuprofen
    }

    smiles = mock_molecules.get(disease_name, "CCO")  # Default: Ethanol
    return jsonify({"molecule": smiles})

@app.route("/predict_efficacy", methods=["POST"])
def predict_efficacy():
    data = request.get_json()
    smiles = data.get("smiles")
    features = get_molecular_features(smiles)

    if features is None:
        return jsonify({"error": "Invalid SMILES"}), 400

    features_df = pd.DataFrame([features])
    prediction = model.predict(features_df)[0]
    return jsonify({"efficacy": round(prediction, 2)})

if __name__ == "__main__":
    app.run(debug=True)

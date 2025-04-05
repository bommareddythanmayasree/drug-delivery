from flask import Flask, request, jsonify, render_template
import joblib
import pandas as pd
import random
import os

app = Flask(__name__, template_folder="templates", static_folder="static")

# Load trained model (replace with your actual model file)

model_path = os.getenv("MODEL_PATH", "drug_discovery_model.joblib")
model = joblib.load(model_path)

# Mock dictionary for molecule generation (using example SMILES strings)
MOCK_MOLECULES = {
    "cancer": "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin (just as a placeholder)
    "diabetes": "CN1C=NC2=C1C(=O)N(C(=O)N2)C",  # Metformin
    "fever": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",   # Ibuprofen
}

# Mock feature extraction (returns random but reasonable dummy values)
def extract_features(smiles):
    if not smiles or not isinstance(smiles, str):
        return None
    return {
        "Molecular_Weight": round(random.uniform(150, 500), 2),
        "LogP": round(random.uniform(-2, 5), 2),
        "H_Bond_Donors": random.randint(0, 5),
        "H_Bond_Acceptors": random.randint(0, 10),
        "TPSA": round(random.uniform(20, 120), 2)
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

    if not smiles or not isinstance(smiles, str):
        return jsonify({"error": "Invalid molecule."}), 400

    # Mock "optimization": reverse the SMILES string just for demonstration
    optimized_smiles = smiles[::-1]
    return jsonify({"optimized": optimized_smiles})

@app.route("/generate_report", methods=["POST"])
def generate_report():
    data = request.get_json()
    smiles = data.get("molecule")
    features = extract_features(smiles)

    if features is None:
        return jsonify({"error": "Invalid molecule."}), 400

    # Generate a simple report as string
    report = "\n".join([f"{k}: {v}" for k, v in features.items()])
    return jsonify({"report": report})

if __name__ == "__main__":
    app.run(debug=True)

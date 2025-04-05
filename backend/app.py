from flask import Flask, request, jsonify, render_template
import joblib
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

# ✅ Load AI Model with Error Handling
model_path = os.path.join(os.getcwd(), "model", "drug_discovery_model.joblib")

if not os.path.exists(model_path):
    raise FileNotFoundError(f"Model file not found at {model_path}. Ensure the model exists in the 'model/' folder.")

model = joblib.load(model_path)
print("✅ Model loaded successfully!")

# ✅ Function to analyze a molecule
def analyze_molecule(smiles):
    if not smiles:
        return {"error": "No SMILES notation provided."}
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES notation."}

    molwt = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)

    features = np.array([[molwt, tpsa, logp, h_donors, h_acceptors]])
    
    try:
        efficacy_pred = model.predict(features)[0]
    except Exception as e:
        return {"error": f"Model prediction failed: {str(e)}"}

    return {
        "Molecular Weight": molwt,
        "TPSA": tpsa,
        "LogP": logp,
        "H Bond Donors": h_donors,
        "H Bond Acceptors": h_acceptors,
        "Predicted Efficacy Score": efficacy_pred
    }

# ✅ Home Route
@app.route("/")
def home():
    return render_template("index.html")

# ✅ Generate Molecule (Placeholder for a real generator)
@app.route("/generate_molecule", methods=["POST"])
def generate_molecule():
    data = request.get_json()
    disease_target = data.get("disease_target", "default_disease")

    # Generate a simple carbon chain molecule (Placeholder for real AI molecule generation)
    mol = Chem.RWMol()
    for _ in range(10):
        mol.AddAtom(Chem.Atom("C"))
    smiles = Chem.MolToSmiles(mol)

    return jsonify({"smiles": smiles, "message": f"Generated molecule for {disease_target}."})

# ✅ Analyze Molecular Properties
@app.route("/analyze", methods=["POST"])
def analyze():
    data = request.get_json()
    smiles = data.get("smiles", "")

    result = analyze_molecule(smiles)
    return jsonify(result)

# ✅ Optimize Molecule Properties
@app.route("/optimize", methods=["POST"])
def optimize():
    data = request.get_json()
    smiles = data.get("smiles", "")

    if not smiles:
        return jsonify({"error": "No SMILES provided."})

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return jsonify({"error": "Invalid molecule."})

    # Example Optimization (Modify this logic)
    molwt = Descriptors.MolWt(mol) * 1.02  # Example change
    logp = Descriptors.MolLogP(mol) * 1.1  # Example change

    try:
        efficacy_pred = model.predict([[molwt, 50, logp, 1, 1]])[0]
    except Exception as e:
        return jsonify({"error": f"Optimization failed: {str(e)}"})

    return jsonify({
        "optimized_smiles": smiles,
        "optimized_eff": efficacy_pred,
        "message": "Optimization completed successfully."
    })

# ✅ Generate Report
@app.route("/report", methods=["POST"])
def report():
    data = request.get_json()
    smiles = data.get("smiles", "")

    if not smiles:
        return jsonify({"error": "No SMILES provided."})

    result = analyze_molecule(smiles)

    if "error" in result:
        return jsonify({"error": result["error"]})

    report_text = f"""
    📜 **Drug Report**
    - **SMILES**: {smiles}
    - **Molecular Weight**: {result["Molecular Weight"]}
    - **TPSA**: {result["TPSA"]}
    - **LogP**: {result["LogP"]}
    - **H Bond Donors**: {result["H Bond Donors"]}
    - **H Bond Acceptors**: {result["H Bond Acceptors"]}
    - **Efficacy Score**: {result["Predicted Efficacy Score"]}
    """

    return jsonify({"report": report_text})

# ✅ Run Flask App
if __name__ == "__main__":
    app.run(debug=True)

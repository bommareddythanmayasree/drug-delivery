from flask import Flask, request, jsonify, render_template
import joblib
import numpy as np
import os
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

#  Load AI Model with Error Handling
model_path = os.path.join(os.getcwd(), "drug_discovery_model.joblib")

if not os.path.exists(model_path):
    raise FileNotFoundError(f"Model file not found at {model_path}. Ensure the model exists in the 'model/' folder.")

model = joblib.load(model_path)
print(" Model loaded successfully!")

#  Function to analyze a molecule
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

#  Home Route
@app.route("/")
def home():
    return render_template("index.html")

@app.route("/generate_molecule", methods=["POST"])
def generate_molecule():
    try:
        if model is None:
            return jsonify({"error": "Model not loaded"}), 500

        data = request.get_json()
        disease_target = data.get("disease_target", "")

        if not disease_target:
            return jsonify({"error": "Disease target is missing"}), 400

        # Predict elements using the model
        elements = model.predict([disease_target])[0]  # e.g., ['C', 'H', 'N', 'O']

        # Join elements to simulate a SMILES string
        smiles = "".join(elements)

        return jsonify({
            "smiles": smiles,
            "message": f"Generated molecule for '{disease_target}'"
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500
#  Analyze Molecular Properties based on disease_target
@app.route("/analyze_molecule", methods=["POST"])
def analyze_molecule():
    data = request.get_json()
    disease_target = data.get("disease_target", "")

    if not disease_target:
        return jsonify({"error": "No disease target provided."})

    try:
        # Step 1: Get elements from model
        elements = model.predict([disease_target])[0]

        # Step 2: Build molecule and convert to SMILES
        mol = Chem.RWMol()
        atom_indices = []

        for element in elements:
            atom_idx = mol.AddAtom(Chem.Atom(element))
            atom_indices.append(atom_idx)

        for i in range(len(atom_indices) - 1):
            mol.AddBond(atom_indices[i], atom_indices[i + 1], Chem.BondType.SINGLE)

        smiles = Chem.MolToSmiles(mol)

        # Step 3: Analyze molecule
        result = analyze_molecule(smiles)
        if "error" in result:
            return jsonify({"error": result["error"]})

        # Step 4: Return results
        return jsonify({
            "smiles": smiles,
            "likelihood": result["Predicted Efficacy Score"],
            "details": result
        })

    except Exception as e:
        return jsonify({"error": f"Failed to analyze molecule: {str(e)}"}), 500



#  Optimize Molecule Properties
@app.route("/optimize_molecule", methods=["POST"])
def optimize_molecule():
    data = request.get_json()
    smiles = data.get("molecule", "")

    if not smiles:
        return jsonify({"error": "No SMILES provided."})

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return jsonify({"error": "Invalid molecule."})

    molwt = Descriptors.MolWt(mol) * 1.02
    logp = Descriptors.MolLogP(mol) * 1.1

    try:
        efficacy_pred = model.predict([[molwt, 50, logp, 1, 1]])[0]
    except Exception as e:
        return jsonify({"error": f"Optimization failed: {str(e)}"})

    return jsonify({
        "optimized": smiles,
        "score": efficacy_pred
    })

#  Generate Report
@app.route("/generate_report", methods=["POST"])
def generate_report():
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

#  Run Flask App
if __name__ == "__main__":
    app.run(debug=True)

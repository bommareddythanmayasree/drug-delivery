from flask import Flask, request, jsonify, render_template
import os
import joblib
import numpy as np
import pandas as pd  # Required if your generate_molecule() returns a DataFrame

app = Flask(__name__)

# Load ML model
try:
    current_directory = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(current_directory, "drug_discovery_model.joblib")
    model = joblib.load(model_path)
except Exception as e:
    model = None
    print(f"Error loading model: {e}")

@app.route("/")
def home():
    return render_template("index.html")

# New: Enhanced Molecule Generator API
@app.route("/generate_molecule", methods=["POST"])
def generate_molecule_api():
    try:
        data = request.get_json()
        disease_target = data.get("disease", "Cancer")  # fallback to 'Cancer' if not provided

        # This function must be defined elsewhere in your code
        generated_smiles, generated_features_df, efficacy_score = generate_molecule(disease_target)

        features_dict = generated_features_df.iloc[0].to_dict()

        return jsonify({
            "smiles": generated_smiles,
            "features": features_dict,
            "efficacy_score": efficacy_score
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/analyze_molecule", methods=["POST"])
def analyze_molecule():
    data = request.get_json()
    molecule = data.get("molecule", "").strip()

    if not molecule:
        return jsonify({"error": "Molecule SMILES is required"}), 400

    if model is None:
        return jsonify({"error": "Model not loaded"}), 500

    try:
        features = [len(molecule)]  # Dummy feature
        prediction = model.predict(np.array(features).reshape(1, -1)).tolist()
        return jsonify({"likelihood": prediction[0]})
    except Exception as e:
        return jsonify({"error": f"Prediction failed: {str(e)}"}), 500

@app.route("/optimize_molecule", methods=["POST"])
def optimize_molecule():
    data = request.get_json()
    molecule = data.get("molecule", "").strip()

    if not molecule:
        return jsonify({"error": "Molecule SMILES is required"}), 400

    optimized = molecule[::-1]
    return jsonify({"optimized": optimized})

@app.route("/generate_report", methods=["POST"])
def generate_report():
    data = request.get_json()
    molecule = data.get("molecule", "").strip()

    if not molecule:
        return jsonify({"error": "Molecule is required"}), 400

    report = (
        f"Report for Molecule:\n"
        f"- SMILES: {molecule}\n"
        f"- Length: {len(molecule)}\n"
        f"- Note: This is a sample report."
    )
    return jsonify({"report": report})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080, debug=True)

from flask import Flask, request, jsonify, render_template
import os
import joblib
import numpy as np
import pandas as pd
import traceback

app = Flask(__name__)

# Load the ML model
model = None
try:
    current_directory = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(current_directory, "drug_discovery_model.joblib")
    model = joblib.load(model_path)
except Exception as e:
    print(f"[ERROR] Failed to load model: {e}")

# Dummy Molecule Generator
def generate_molecule(disease_target):
    smiles = f"SMILES_{disease_target}"
    features_df = pd.DataFrame([{
        "length": len(smiles),
        "dummy_feature": 1.0
    }])
    efficacy = 0.87
    return smiles, features_df, efficacy


@app.route("/")
def home():
    return render_template("index.html")


@app.route("/generate_molecule", methods=["POST"])
def generate_molecule_api():
    try:
        data = request.get_json()
        disease_target = data.get("disease", "Cancer")

        smiles, features_df, efficacy_score = generate_molecule(disease_target)
        features_dict = features_df.iloc[0].to_dict()

        return jsonify({
            "smiles": smiles,
            "features": features_dict,
            "efficacy_score": efficacy_score
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/analyze_molecule", methods=["POST"])
def analyze_molecule():
    try:
        data = request.get_json()
        molecule = data.get("molecule", "").strip()

        if not molecule:
            return jsonify({"error": "Molecule SMILES is required"}), 400
        if model is None:
            return jsonify({"error": "Model not loaded"}), 500

        features_df = pd.DataFrame([{"length": len(molecule)}])
        prediction = model.predict(features_df).tolist()

        return jsonify({"likelihood": prediction[0]})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Prediction failed: {str(e)}"}), 500


@app.route("/optimize_molecule", methods=["POST"])
def optimize_molecule():
    try:
        data = request.get_json()
        molecule = data.get("molecule", "").strip()

        if not molecule:
            return jsonify({"error": "Molecule is required"}), 400

        optimized = molecule[::-1]  # Dummy logic: reverse the SMILES string
        return jsonify({"optimized": optimized})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/generate_report", methods=["POST"])
def generate_report():
    try:
        data = request.get_json()
        molecule = data.get("molecule", "").strip()

        if not molecule:
            return jsonify({"error": "Molecule is required"}), 400

        report = f"Report for Molecule:\n- SMILES: {molecule}\n- Length: {len(molecule)}"
        return jsonify({"report": report})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# Do not use Flask dev server in production; gunicorn is used instead.
if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)

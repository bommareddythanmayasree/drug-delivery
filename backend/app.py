from flask import Flask, request, jsonify, render_template
import os
import joblib
import numpy as np

app = Flask(__name__)

# Load your ML model
current_directory = os.path.dirname(os.path.abspath(__file__))
model_path = os.path.join(current_directory, "drug_discovery_model.joblib")
model = joblib.load(model_path)

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/generate_molecule", methods=["POST"])
def generate_molecule():
    data = request.get_json()
    disease_target = data.get("disease_target", "").lower()

    # Dummy molecule generation logic – replace with your actual method
    molecule = f"Generated_SMILES_for_{disease_target}"
    return jsonify({"molecule": molecule})

@app.route("/analyze_molecule", methods=["POST"])
def analyze_molecule():
    data = request.get_json()
    molecule = data.get("molecule", "")

    # Dummy analysis logic – replace with actual analysis
    # For demo: pretend the model uses simple numeric features (e.g., length of molecule string)
    features = [len(molecule)]
    prediction = model.predict(np.array(features).reshape(1, -1)).tolist()
    
    return jsonify({"likelihood": prediction[0]})

@app.route("/optimize_molecule", methods=["POST"])
def optimize_molecule():
    data = request.get_json()
    molecule = data.get("molecule", "")

    # Dummy optimization logic – replace with real algorithm
    optimized = molecule[::-1]  # Just reverse string as a placeholder
    return jsonify({"optimized": optimized})

@app.route("/generate_report", methods=["POST"])
def generate_report():
    data = request.get_json()
    molecule = data.get("molecule", "")

    # Dummy report generation
    report = f"Report for Molecule:\n- SMILES: {molecule}\n- Length: {len(molecule)}"
    return jsonify({"report": report})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)

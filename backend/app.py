from flask import Flask, request, jsonify, render_template
import os
import joblib
import numpy as np

app = Flask(__name__)  # Initialize the Flask app

# Get the absolute path of the current script
current_directory = os.path.dirname(os.path.abspath(__file__))

# Construct the correct path to the model file
model_path = os.path.join(current_directory, "drug_discovery_model.joblib")

# Load the model
model = joblib.load(model_path)

# Home route
@app.route("/")
def home():
    return render_template("index.html")

# Prediction API
@app.route("/predict", methods=["POST"])
def predict():
    data = request.get_json()
    features = np.array(data["features"]).reshape(1, -1)
    prediction = model.predict(features).tolist()
    return jsonify({"prediction": prediction})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)





import pandas as pd
import numpy as np
import joblib
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
import matplotlib.pyplot as plt
from tqdm import tqdm
import random

# 📌 Step 1: Load & Preprocess the Dataset
file_path = "/content/Updated_Drug_Discovery_Dataset.csv"
df = pd.read_csv(file_path)
df = df.dropna()

# Select relevant features
required_columns = ["Molecular_Weight", "TPSA", "LogP", "H_Bond_Donors", "H_Bond_Acceptors", "Efficacy_Score"]
df = df[required_columns]

# Define input (X) and output (y) variables
X = df.drop(columns=["Efficacy_Score"])
y = df["Efficacy_Score"]

# 📌 Step 2: Split Data into Training & Testing Sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 📌 Step 3: Train the AI Model (Random Forest)
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Save the model
model_filename = "drug_discovery_model.joblib"
joblib.dump(model, model_filename)

# 📌 Step 4: Molecule Generation
def generate_molecule(disease_target):
    print(f"Generating molecule for disease: {disease_target}")

    # Generate random molecular features
    mol_features = np.random.uniform(low=1, high=10, size=(1, 5))
    feature_names = ["Molecular_Weight", "TPSA", "LogP", "H_Bond_Donors", "H_Bond_Acceptors"]
    mol_features_df = pd.DataFrame(mol_features, columns=feature_names)

    # Predict efficacy score
    predicted_eff = model.predict(mol_features_df)[0]

    # Generate a molecule using RDKit
    mol = Chem.RWMol()
    atoms = ["C", "N", "O", "S", "Cl", "Br"]  # Common atoms in drugs
    for _ in range(random.randint(5, 15)):  # Generate random molecule length
        atom = random.choice(atoms)
        mol.AddAtom(Chem.Atom(atom))

    # Convert to SMILES
    smiles = Chem.MolToSmiles(mol)

    # Ensure the molecule is valid
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        smiles = "Invalid molecule generated, try again."

    return smiles, mol_features_df, predicted_eff

# 📌 Step 5: Molecular Analysis & Efficacy Prediction
def analyze_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES notation"

    molwt = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)

    features = np.array([[molwt, tpsa, logp, h_donors, h_acceptors]])
    efficacy_pred = model.predict(features)[0]

    return {
        "Molecular Weight": molwt,
        "TPSA": tpsa,
        "LogP": logp,
        "H Bond Donors": h_donors,
        "H Bond Acceptors": h_acceptors,
        "Predicted Efficacy Score": efficacy_pred
    }

def optimize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Invalid SMILES notation")
        return None, None, None

    # Extract molecular properties
    molwt = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)

    # Store original features
    original_features = np.array([[molwt, tpsa, logp, h_donors, h_acceptors]])
    original_efficacy = model.predict(original_features)[0]  # Predict original efficacy

    # 🔹 Smart Optimization Strategy:
    optimized_features = original_features.copy()
    optimized_features[0, 2] *= 1.1  # Increase LogP by 10%
    optimized_features[0, 1] *= 1.05  # Increase TPSA slightly
    optimized_features[0, 0] *= 1.02  # Small increase in Molecular Weight
    optimized_features[0, 3] = max(0, optimized_features[0, 3] - 1)  # Reduce H-Bond Donors
    optimized_features[0, 4] = max(0, optimized_features[0, 4] - 1)  # Reduce H-Bond Acceptors

    # Predict new efficacy after optimization
    optimized_efficacy = model.predict(optimized_features)[0]

    # If optimization made efficacy worse, revert to the original
    if optimized_efficacy < original_efficacy:
        optimized_features = original_features.copy()
        optimized_efficacy = original_efficacy

    optimized_smiles = Chem.MolToSmiles(mol)  # Keeping original SMILES

    print("\n🔹 **Optimization & Validation Report** 🔹")
    print(f"Original Molecule: {smiles}")
    print(f"Optimized Molecule: {optimized_smiles}")
    print("\n📌 **Molecular Properties Before & After Optimization** 📌")
    print(f"- Molecular Weight: {molwt:.2f} ➝ {optimized_features[0, 0]:.2f}")
    print(f"- TPSA: {tpsa:.2f} ➝ {optimized_features[0, 1]:.2f}")
    print(f"- LogP: {logp:.2f} ➝ {optimized_features[0, 2]:.2f}")
    print(f"- H-Bond Donors: {h_donors} ➝ {int(optimized_features[0, 3])}")
    print(f"- H-Bond Acceptors: {h_acceptors} ➝ {int(optimized_features[0, 4])}")
    print(f"\n✅ **Efficacy Score Before Optimization:** {original_efficacy:.4f}")
    print(f"✅ **Efficacy Score After Optimization:** {optimized_efficacy:.4f}")

    return optimized_smiles, optimized_features, optimized_efficacy

def calculate_solubility(smiles):
    """Estimate solubility based on molecular properties."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES"

    logp = Descriptors.MolLogP(mol)
    molwt = Descriptors.MolWt(mol)
    logS = -0.75 * logp - 0.006 * molwt + 0.1

    if logS > -2:
        return f"High (LogS: {logS:.2f})"
    elif logS > -4:
        return f"Medium (LogS: {logS:.2f})"
    else:
        return f"Low (LogS: {logS:.2f})"

def calculate_toxicity(smiles):
    """Simple toxicity prediction based on molecular properties."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES"

    logp = Descriptors.MolLogP(mol)
    molwt = Descriptors.MolWt(mol)

    if logp > 5 or molwt > 500:
        return "High Toxicity"
    elif logp > 3 or molwt > 350:
        return "Medium Toxicity"
    else:
        return "Low Toxicity"

def calculate_druglikeness(smiles):
    """Evaluate drug-likeness based on Lipinski's Rule of Five."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES"

    # Get molecular properties
    molwt = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)

    # Lipinski's Rule of Five criteria
    violations = 0
    if molwt > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if h_donors > 5:
        violations += 1
    if h_acceptors > 10:
        violations += 1

    # Additional criteria (Veber's rules)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)

    if rotatable_bonds > 10:
        violations += 1
    if tpsa > 140:
        violations += 1

    # Determine drug-likeness
    if violations == 0:
        return "Excellent,Score -> 5/5 (0 violations)"
    elif violations == 1:
        return "Good,Score -> 4/5 (1 violation)"
    elif violations == 2:
        return "Moderate,Score -> 3/5 (2 violations)"
    else:
        return "Poor,Score less than 2 (3+ violations)"

def generate_report(smiles):
    """Generate a full drug discovery report with all analyses."""
    analysis = analyze_molecule(smiles)

    if isinstance(analysis, str):  # Check for invalid SMILES
        return analysis

    # Calculate additional properties
    toxicity = calculate_toxicity(smiles)
    solubility = calculate_solubility(smiles)
    druglikeness = calculate_druglikeness(smiles)

    # Generate report
    report = f"""
    📜 **Drug Discovery Report** 📜
    --------------------------------
    🧪 **Molecular Analysis**
    - **SMILES Representation:** {smiles}
    - **Molecular Weight:** {analysis["Molecular Weight"]:.2f}
    - **TPSA (Topological Polar Surface Area):** {analysis["TPSA"]:.2f}
    - **LogP (Lipophilicity):** {analysis["LogP"]:.2f}
    - **H-Bond Donors:** {analysis["H Bond Donors"]}
    - **H-Bond Acceptors:** {analysis["H Bond Acceptors"]}

    🔥 **Predicted Efficacy Score:** {analysis["Predicted Efficacy Score"]:.4f}
    ☠ **Estimated Toxicity:** {toxicity}
    💧 **Estimated Solubility:** {solubility}
    💊 **Drug-likeness (Lipinski's Rule of Five):** {druglikeness}
    --------------------------------
    """

    return report
joblib.dump(model, "drug_discovery_model.joblib")


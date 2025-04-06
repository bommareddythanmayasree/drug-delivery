document.addEventListener("DOMContentLoaded", function () {
    const generateBtn = document.getElementById("generateMolecule");
    const analyzeBtn = document.getElementById("analyzeMolecule");
    const optimizeBtn = document.getElementById("optimizeMolecule");
    const reportBtn = document.getElementById("generateReport");

    const resultDiv = document.getElementById("result");

    // API Base URL
    const API_BASE_URL = "http://127.0.0.1:5000";  // Ensure Flask is running locally

    // Generate Molecule
    generateBtn.addEventListener("click", async function () {
        const diseaseTarget = document.getElementById("diseaseInput").value;
        if (!diseaseTarget) {
            alert("Please enter a disease target!");
            return;
        }

        const response = await fetch(`${API_BASE_URL}/generate_molecule`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ disease_target: diseaseTarget })
        });

        const data = await response.json();
        resultDiv.innerHTML = `<h3>Generated Molecule:</h3><p>${data.molecule}</p>`;
    });

    // Molecular Analysis
    analyzeBtn.addEventListener("click", async function () {
        const moleculeStructure = document.getElementById("moleculeInput").value;
        if (!moleculeStructure) {
            alert("Please enter a molecule structure!");
            return;
        }

        const response = await fetch(`${API_BASE_URL}/analyze_molecule`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ molecule: moleculeStructure })
        });

        const data = await response.json();
        resultDiv.innerHTML = `<h3>Analysis Results:</h3>
                               <p>Drug Likelihood: ${data.drug_likeliness}</p>
                               <p>Efficacy Score: ${data.efficacy_score}</p>`;
    });

    // Optimize Molecule
    optimizeBtn.addEventListener("click", async function () {
        const moleculeStructure = document.getElementById("moleculeInput").value;
        if (!moleculeStructure) {
            alert("Please enter a molecule structure!");
            return;
        }

        const response = await fetch(`${API_BASE_URL}/optimize_molecule`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ molecule: moleculeStructure })
        });

        const data = await response.json();
        resultDiv.innerHTML = `<h3>Optimized Molecule:</h3><p>${data.optimized_molecule}</p>`;
    });

    // Generate Report
    reportBtn.addEventListener("click", async function () {
        const response = await fetch(`${API_BASE_URL}/generate_report`, {
            method: "GET",
            headers: { "Content-Type": "application/json" }
        });

        const data = await response.json();
        resultDiv.innerHTML = `<h3>Report:</h3><pre>${data.report}</pre>`;
    });
});

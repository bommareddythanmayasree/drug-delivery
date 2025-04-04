document.addEventListener("DOMContentLoaded", function () {
    const generateBtn = document.getElementById("generateMolecule");
    const analyzeBtn = document.getElementById("analyzeMolecule");
    const optimizeBtn = document.getElementById("optimizeMolecule");
    const reportBtn = document.getElementById("generateReport");

    const resultDiv = document.getElementById("result");

    // Change this when deploying
    const API_BASE_URL = "https://your-backend.onrender.com";  

    // Helper function for API calls
    async function fetchAPI(endpoint, method, body = null) {
        try {
            const response = await fetch(`${API_BASE_URL}/${endpoint}`, {
                method: method,
                headers: { "Content-Type": "application/json" },
                body: body ? JSON.stringify(body) : null,
            });

            if (!response.ok) {
                throw new Error(`Error: ${response.status} ${response.statusText}`);
            }

            return await response.json();
        } catch (error) {
            alert("Failed to connect to the server. Please try again later.");
            console.error("API Error:", error);
        }
    }

    // Generate Molecule
    generateBtn.addEventListener("click", async function () {
        const diseaseTarget = document.getElementById("diseaseInput").value.trim();
        if (!diseaseTarget) {
            alert("Please enter a disease target!");
            return;
        }

        const data = await fetchAPI("generate_molecule", "POST", { disease_target: diseaseTarget });
        if (data) {
            resultDiv.innerHTML = `<h3>Generated Molecule:</h3><p>${data.molecule}</p>`;
        }
    });

    // Molecular Analysis
    analyzeBtn.addEventListener("click", async function () {
        const moleculeStructure = document.getElementById("moleculeInput").value.trim();
        if (!moleculeStructure) {
            alert("Please enter a molecule structure!");
            return;
        }

        const data = await fetchAPI("analyze_molecule", "POST", { molecule: moleculeStructure });
        if (data) {
            resultDiv.innerHTML = `<h3>Analysis Results:</h3>
                                   <p><strong>Drug Likelihood:</strong> ${data.drug_likeliness}</p>
                                   <p><strong>Efficacy Score:</strong> ${data.efficacy_score}</p>`;
        }
    });

    // Optimize Molecule
    optimizeBtn.addEventListener("click", async function () {
        const moleculeStructure = document.getElementById("moleculeInput").value.trim();
        if (!moleculeStructure) {
            alert("Please enter a molecule structure!");
            return;
        }

        const data = await fetchAPI("optimize_molecule", "POST", { molecule: moleculeStructure });
        if (data) {
            resultDiv.innerHTML = `<h3>Optimized Molecule:</h3><p>${data.optimized_molecule}</p>`;
        }
    });

    // Generate Report
    reportBtn.addEventListener("click", async function () {
        const data = await fetchAPI("generate_report", "GET");
        if (data) {
            resultDiv.innerHTML = `<h3>Report:</h3><pre>${data.report}</pre>`;
        }
    });
});

document.addEventListener("DOMContentLoaded", function () {
    const generateBtn = document.getElementById("generateMolecule");
    const analyzeBtn = document.getElementById("analyzeMolecule");
    const optimizeBtn = document.getElementById("optimizeMolecule");
    const reportBtn = document.getElementById("generateReport");

    const resultDiv = document.getElementById("result");

    //  Set correct backend URL
    const API_BASE_URL = "https://drug-delivery.onrender.com";  

    //  Helper function for API calls with improved error handling
    async function fetchAPI(endpoint, method, body = null) {
        try {
            const response = await fetch(`${API_BASE_URL}/${endpoint}`, {
                method: method,
                headers: { "Content-Type": "application/json" },
                body: body ? JSON.stringify(body) : null,
            });

            if (!response.ok) {
                throw new Error(`Server Error: ${response.status} ${response.statusText}`);
            }

            return await response.json();
        } catch (error) {
            resultDiv.innerHTML = `<p style="color: red;"><strong>Error:</strong> ${error.message}</p>`;
            console.error("API Error:", error);
        }
    }

    //  Generate Molecule
    generateBtn.addEventListener("click", async function () {
        const diseaseTarget = document.getElementById("diseaseInput").value.trim();
        if (!diseaseTarget) {
            resultDiv.innerHTML = `<p style="color: red;">Please enter a disease target!</p>`;
            return;
        }

        const data = await fetchAPI("generate_molecule", "POST", { disease_target: diseaseTarget });
        if (data) {
            resultDiv.innerHTML = `<h3>Generated Molecule:</h3><p>${data.molecule}</p>`;
        }
    });

    //  Molecular Analysis
    analyzeBtn.addEventListener("click", async function () {
        const moleculeStructure = document.getElementById("moleculeInput").value.trim();
        if (!moleculeStructure) {
            resultDiv.innerHTML = `<p style="color: red;">Please enter a molecule structure!</p>`;
            return;
        }

        const data = await fetchAPI("analyze_molecule", "POST", { molecule: moleculeStructure });
        if (data) {
            resultDiv.innerHTML = `<h3>Analysis Results:</h3>
                                   <p><strong>Drug Likelihood:</strong>

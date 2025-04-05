document.addEventListener("DOMContentLoaded", function () {
    const API_BASE_URL = "https://drug-delivery.onrender.com";  // ✅ Replace with your actual backend URL

    // Button elements
    const generateBtn = document.querySelector("#generate button");
    const analyzeBtn = document.querySelector("#analyze button");
    const optimizeBtn = document.querySelector("#optimize button");
    const reportBtn = document.querySelector("#report button");

    // Result areas
    const generatedOutput = document.getElementById("generated_molecule");
    const analysisOutput = document.getElementById("analysis_results");
    const optimizedOutput = document.getElementById("optimized_molecule");
    const reportOutput = document.getElementById("report_message");

    // Tab switching
    window.showTab = function (tabId) {
        document.querySelectorAll(".tab-content").forEach(tab => {
            tab.style.display = "none";
        });
        document.getElementById(tabId).style.display = "block";

        document.querySelectorAll(".tab-btn").forEach(btn => {
            btn.classList.remove("active");
        });
        document.querySelector(`[onclick="showTab('${tabId}')"]`).classList.add("active");
    };

    // Helper for API calls
    async function fetchAPI(endpoint, method, body = null) {
        try {
            const response = await fetch(`${API_BASE_URL}/${endpoint}`, {
                method: method,
                headers: { "Content-Type": "application/json" },
                body: body ? JSON.stringify(body) : null
            });

            if (!response.ok) {
                throw new Error(`Server Error: ${response.status} ${response.statusText}`);
            }

            return await response.json();
        } catch (error) {
            return { error: error.message };
        }
    }

    // 🔬 Generate Molecule
    generateBtn.addEventListener("click", async function () {
        const target = document.getElementById("disease_target").value.trim();
        if (!target) {
            generatedOutput.innerHTML = `<span style="color:red;">Please enter a disease target.</span>`;
            return;
        }

        const data = await fetchAPI("generate_molecule", "POST", { disease_target: target });
        if (data.error) {
            generatedOutput.innerHTML = `<span style="color:red;">Error: ${data.error}</span>`;
        } else {
            generatedOutput.innerHTML = `Generated Molecule: <strong>${data.molecule}</strong>`;
        }
    });

    // 🧪 Molecular Analysis
    analyzeBtn.addEventListener("click", async function () {
        const smiles = document.getElementById("molecule_input").value.trim();
        if (!smiles) {
            analysisOutput.innerHTML = `<span style="color:red;">Please enter a molecule SMILES string.</span>`;
            return;
        }

        const data = await fetchAPI("analyze_molecule", "POST", { molecule: smiles });
        if (data.error) {
            analysisOutput.innerHTML = `<span style="color:red;">Error: ${data.error}</span>`;
        } else {
            analysisOutput.innerHTML = `Drug Likelihood: <strong>${data.likelihood}</strong>`;
        }
    });

    // 🛠 Optimize Molecule
    optimizeBtn.addEventListener("click", async function () {
        const input = document.getElementById("optimize_input").value.trim();
        if (!input) {
            optimizedOutput.innerHTML = `<span style="color:red;">Please enter a molecule SMILES string.</span>`;
            return;
        }

        const data = await fetchAPI("optimize_molecule", "POST", { molecule: input });
        if (data.error) {
            optimizedOutput.innerHTML = `<span style="color:red;">Error: ${data.error}</span>`;
        } else {
            optimizedOutput.innerHTML = `Optimized Molecule: <strong>${data.optimized}</strong>`;
        }
    });

    // 📜 Generate Report
    reportBtn.addEventListener("click", async function () {
        const input = document.getElementById("report_molecule").value.trim();
        if (!input) {
            reportOutput.innerHTML = `<span style="color:red;">Please enter a molecule to generate report.</span>`;
            return;
        }

        const data = await fetchAPI("generate_report", "POST", { molecule: input });
        if (data.error) {
            reportOutput.innerHTML = `<span style="color:red;">Error: ${data.error}</span>`;
        } else {
            reportOutput.innerHTML = `Report Generated: <pre>${data.report}</pre>`;
        }
    });

    // Load default tab
    showTab("generate");
});

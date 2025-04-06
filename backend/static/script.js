document.addEventListener("DOMContentLoaded", function () {
    const API_BASE_URL = "https://drug-delivery.onrender.com";

    // Tab logic
    const tabs = document.querySelectorAll(".tab-btn");
    tabs.forEach(tab => {
        tab.addEventListener("click", () => {
            document.querySelectorAll(".tab-content").forEach(tab => tab.style.display = "none");
            document.getElementById(tab.dataset.tab).style.display = "block";
            tabs.forEach(btn => btn.classList.remove("active"));
            tab.classList.add("active");
        });
    });

    // Helper for API calls
    async function fetchAPI(endpoint, method, body = null) {
        try {
            console.log("Sending request to:", `${API_BASE_URL}/${endpoint}`);
            const response = await fetch(`${API_BASE_URL}/${endpoint}`, {
                method,
                mode: 'cors', // Ensures cross-origin requests work properly
                headers: {
                    "Content-Type": "application/json"
                },
                body: body ? JSON.stringify(body) : null
            });

            if (!response.ok) throw new Error(`Server Error: ${response.status} ${response.statusText}`);
            return await response.json();
        } catch (error) {
            console.error("Fetch error:", error); // Debug log
            return { error: error.message };
        }
    }

    // GENERATE MOLECULE
    document.querySelector("#generate button").addEventListener("click", async () => {
        const target = document.getElementById("disease_target").value.trim();
        const output = document.getElementById("generated_molecule");

        if (!target) {
            output.innerHTML = `<span style="color:red;">Please enter a disease target.</span>`;
            return;
        }

        const data = await fetchAPI("generate_molecule", "POST", { disease_target: target });
        output.innerHTML = data.error
            ? `<span style="color:red;">Error: ${data.error}</span>`
            : `Generated Molecule: <strong>${data.smiles}</strong>`;
    });

    // ANALYZE MOLECULE
    document.querySelector("#analyze button").addEventListener("click", async () => {
        const smiles = document.getElementById("molecule_input").value.trim();
        const output = document.getElementById("analysis_results");

        if (!smiles) {
            output.innerHTML = `<span style="color:red;">Please enter a molecule SMILES string.</span>`;
            return;
        }

        const data = await fetchAPI("analyze_molecule", "POST", { disease_target: smiles });
        output.innerHTML = data.error
            ? `<span style="color:red;">Error: ${data.error}</span>`
            : `
                <strong>SMILES:</strong> ${data.smiles}<br>
                <strong>Efficacy:</strong> ${data.likelihood}<br>
                <pre>${JSON.stringify(data.details, null, 2)}</pre>
            `;
    });

    // OPTIMIZE MOLECULE
    document.querySelector("#optimize button").addEventListener("click", async () => {
        const input = document.getElementById("optimize_input").value.trim();
        const output = document.getElementById("optimized_molecule");

        if (!input) {
            output.innerHTML = `<span style="color:red;">Please enter a molecule SMILES string.</span>`;
            return;
        }

        const data = await fetchAPI("optimize_molecule", "POST", { molecule: input });
        output.innerHTML = data.error
            ? `<span style="color:red;">Error: ${data.error}</span>`
            : `Optimized Molecule: <strong>${data.optimized}</strong> <br>Score: ${data.score}`;
    });

    // GENERATE REPORT
    document.querySelector("#report button").addEventListener("click", async () => {
        const input = document.getElementById("report_molecule").value.trim();
        const output = document.getElementById("report_message");

        if (!input) {
            output.innerHTML = `<span style="color:red;">Please enter a molecule to generate report.</span>`;
            return

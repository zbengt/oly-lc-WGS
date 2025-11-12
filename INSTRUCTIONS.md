# 🧠 AGENT INSTRUCTIONS

These instructions define how an automated agent (e.g., ChatGPT, Claude) should work within this repository.

---

## 📁 Repository Structure

```
repo/
├── code/      → all analysis or processing scripts
├── data/      → raw and reference input data (read-only)
└── output/    → generated results organized by script
```

**Key convention:**  
Each file in `code/` is prefixed with a sequential number (e.g., `01_`, `02_`, `03_`), indicating **execution order**.  
Each script writes outputs into a **matching subdirectory** inside `output/`.

| Code File | Output Directory | Purpose |
|------------|------------------|----------|
| `code/01_load_data.py` | `output/01_load_data/` | Load and clean input data |
| `code/02_preprocess.R` | `output/02_preprocess/` | Transform and filter data |
| `code/03_analyze.ipynb` | `output/03_analyze/` | Model and visualize results |

---

## ⚙️ Execution Rules

1. **Run code files in ascending numeric order.**  
   Example: run `01_` before `02_`, `02_` before `03_`.

2. **Each step depends on the previous one.**  
   Inputs for `02_` come from `output/01_*/`, etc.

3. **Each script writes results to a matching output subdirectory.**  
   Example:  
   ```python
   os.makedirs("output/02_preprocess", exist_ok=True)
   ```

4. **Never modify files in `data/`.**  
   Treat this directory as read-only input.

---

## 🔁 Input/Output Conventions

### Inputs
- `data/` directory for raw or reference data  
- `output/<previous_step>/` for intermediate results  

### Outputs
- Write all results, plots, and logs to:
  ```
  output/<current_step_name>/
  ```

### File Examples
```
output/02_preprocess/processed_data.csv
output/02_preprocess/metadata.json
output/03_analyze/figures/accuracy_plot.png
```

---

## 🧩 Metadata and Logging

Each script should generate a `metadata.json` file inside its output directory.

**Example:**
```json
{
  "script": "02_preprocess.R",
  "input": "../output/01_load_data/",
  "parameters": {"filter_threshold": 0.05},
  "date": "2025-11-12T12:00:00Z",
  "runtime_seconds": 128.4,
  "random_seed": 42
}
```

Recommended contents:
- script name  
- input sources  
- parameters or arguments  
- execution date/time  
- runtime and software versions  
- notes or comments (optional)

---

## 🧬 Reproducibility Rules

- Each step must be runnable independently once previous outputs exist.  
- Use **relative paths only** (no absolute paths).  
- Always set **random seeds** for reproducibility.  
- Never overwrite outputs unless explicitly rerunning a step.  
- Log parameters and environment info for each execution.  

---

## 🪶 Agent Behavior Protocol

When given an instruction, the agent should:

1. Identify the correct script in `code/` based on task description or prefix.  
2. Verify required inputs exist in `data/` or the previous `output/` directory.  
3. Create the corresponding output subdirectory if needed.  
4. Execute or simulate the logic of the script.  
5. Write all generated files to `output/<script_name>/`.  
6. Save or update the `metadata.json`.  
7. Avoid modifying any data or outputs not belonging to the current step.

---

## 🧱 Naming and Structure

- **Scripts:** `NN_description.ext` (e.g., `03_visualize_results.R`)  
- **Output folders:** `output/NN_description/`  
- **Output files:** use clear, lowercase, hyphenated names (e.g., `summary-stats.csv`, not `SummaryStats.csv`).

Avoid spaces or uppercase letters in file and folder names.

---

## 🚨 Error Handling

If an expected input file or folder is missing:
1. Report which dependency is missing (e.g., `output/01_load_data/clean_data.csv not found`).
2. Do not create placeholder outputs.
3. Stop execution and notify that the previous step must be completed first.

If an output directory already exists:
- Write new files with version suffixes (e.g., `_v2`, `_rerun`), or
- Append a timestamp to avoid overwriting existing results.

---

## 📊 Example Directory Snapshot

```
repo/
├── code/
│   ├── 01_load_data.py
│   ├── 02_preprocess.R
│   ├── 03_visualize_results.R
│   └── 04_analyze.ipynb
├── data/
│   ├── raw_input.csv
│   └── metadata.yaml
└── output/
    ├── 01_load_data/
    │   ├── clean_data.csv
    │   └── metadata.json
    ├── 02_preprocess/
    │   ├── processed_data.csv
    │   └── metadata.json
    ├── 03_visualize_results/
    │   ├── accuracy_plot.png
    │   └── metadata.json
    └── 04_analyze/
        ├── model_output.pkl
        └── summary_report.html
```

---

## ✅ Quick Checklist

- [ ] Execute scripts in numeric order  
- [ ] Read only from `data/` or previous outputs  
- [ ] Write to `output/<script_name>/`  
- [ ] Include a `metadata.json` in every output directory  
- [ ] Avoid overwriting or altering data files  
- [ ] Maintain relative paths and reproducibility

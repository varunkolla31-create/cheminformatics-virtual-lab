# Cheminformatics Virtual Laboratory 🧪

An interactive, modern, beginner-friendly **Streamlit multi-page web application** designed for undergraduate chemistry students to learn fundamental concepts in cheminformatics.

## Overview & Learning Objectives

Students exploring the **Cheminformatics Virtual Laboratory** will master:
- **SMILES notation**: Writing and interpreting Simplified Molecular Input Line Entry System strings.
- **Molecular representation**: Understanding how chemical structures map to computational graphs.
- **2D visualization**: Generating publication-quality 2D chemical depictions with RDKit.
- **3D visualization**: Building, protonating, and energy-minimizing 3D conformers with interactive `py3Dmol` / `stmol` rendering.
- **Molecular descriptors**: Calculating physicochemical descriptors (MW, LogP, TPSA, HBD, HBA, rotatable bonds, rings) and evaluating Lipinski's Rule of Five.
- **Structure-property relationships (SAR/QSPR)**: Comparing multiple compounds via interactive Plotly bar and scatter charts.

## Project Structure

| File | Description |
| :--- | :--- |
| `app.py` | Main Streamlit application with modular chemistry helpers, sidebar navigation (`streamlit-option-menu`), and 7 interactive lab pages |
| `requirements.txt` | Streamlit Community Cloud compatible Python dependencies (uses `rdkit`, not `rdkit-pypi`) |
| `.streamlit/config.toml` | Custom academic theme configuration (`#1E3A8A` primary navy, `#F8FAFC` slate background, `#FFFFFF` cards) |
| `README.md` | Project documentation, local setup guide, and Streamlit Cloud deployment instructions |

## Application Pages

1. **Home**: Hero banner with "Start Learning" quick-action, interactive feature grid, reference table of 6 benchmark molecules, and academic footer.
2. **Theory**: Expandable modules covering SMILES syntax, Molecular Graphs, Molecular Descriptors, and Lipinski's Rule of Five with live examples.
3. **Molecular Visualization**: Dual-tab (`2D | 3D`) viewer with one-click example buttons, friendly SMILES error diagnostics, RDKit 2D drawings, and MMFF94-optimized 3D models (`stick`, `sphere`, `line` styles and customizable color schemes).
4. **Molecular Properties Report**: Physicochemical descriptor calculation displayed as `st.metric` cards, Lipinski pass/fail badge, structured table, and downloadable **CSV** and **PDF** reports stamped with student metadata.
5. **Structure-Property Analysis**: Batch multi-SMILES comparison table, customizable descriptor bar chart, and interactive **MW vs. LogP** chemical space scatter plot.
6. **Assessment**: 10-question interactive multiple-choice quiz (one question at a time) with progress bar, instant feedback explanations, and final achievement badge.
7. **Final Project Report**: Consolidated lab portfolio summarizing student registration info, all molecules studied during the session, descriptor summary table, assessment score, and one-click **PDF** portfolio export.

## Local Installation & Usage

1. **Clone the repository**:
   ```bash
   git clone https://github.com/varunkolla31-create/cheminformatics-virtual-lab.git
   cd cheminformatics-virtual-lab
   ```

2. **Create a virtual environment (Python 3.11 recommended)**:
   ```bash
   python3.11 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Run the Streamlit app**:
   ```bash
   streamlit run app.py
   ```

## Deploying to Streamlit Community Cloud

1. Push `app.py`, `requirements.txt`, `.streamlit/config.toml`, and `README.md` to a public or private GitHub repository.
2. Sign in to [share.streamlit.io](https://share.streamlit.io) and click **New app**.
3. Select your repository (`varunkolla31-create/cheminformatics-virtual-lab`), branch (`main`), and set **Main file path** to `app.py`.
4. Under **Advanced settings**, select **Python 3.11** and click **Deploy!**.

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

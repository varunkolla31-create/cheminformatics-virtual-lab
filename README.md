# 🧪 Cheminformatics Virtual Laboratory (Apple UI/UX Edition)

An Apple.com-inspired minimalist **Streamlit multipage web application** for computational chemistry, molecular visualisation, ADMET property calculation, Tanimoto fingerprint similarity screening, and chemical space exploration.

**Author:** Varun K (`RA2611026050122` | Section: `AIML B`)  
**Institution:** SRM Institute of Science and Technology

---

## ✨ Key Features & Multipage Architecture

| File / Page | Description |
| :--- | :--- |
| `app.py` | **Home:** Apple-style hero banner, quick-navigation pill links, 6 capability cards, and a 10-compound pharmaceutical reference table. |
| `pages/1_Molecule_Viewer.py` | **Molecule Viewer:** SMILES validation, RDKit 2D vector structure depiction, and interactive MMFF94 force-field optimised 3D conformers (`py3Dmol` / `stmol`) with PNG & SDF downloads. |
| `pages/2_Property_Calculator.py` | **Property Calculator:** Computes MW, Crippen LogP, TPSA, HBD, HBA, Rotatable Bonds, Rings, and Fraction Csp3 with a Lipinski Rule-of-Five radar chart and CSV/PDF report exports. |
| `pages/3_Similarity_Search.py` | **Similarity Search:** Screens a 24-compound drug library using Morgan (ECFP4), MACCS Keys (166-bit), and RDKit topological fingerprints ranked by Tanimoto similarity. |
| `pages/4_Dataset_Explorer.py` | **Dataset Explorer:** Supports custom CSV and SDF file uploads, interactive MW/LogP/Lipinski sliders, and Plotly chemical space scatter/bar visualisations. |
| `pages/5_About.py` | **Theory, Quiz & Portfolio:** Expandable cheminformatics theory modules, a 10-question interactive self-assessment quiz, and an official ReportLab PDF portfolio generator. |
| `utils.py` | **Core Engine & Design System:** Cached RDKit/py3Dmol helpers, Apple.com CSS injection (`#FBFBFD` background, `#FFFFFF` cards, `#0071E3` accent, `#D2D2D7` borders), sticky glassmorphic top nav, and mandatory footer. |

---

## 🎨 Apple-Style UI/UX Specifications

- **Palette:** Background `#FBFBFD`, Surface `#FFFFFF`, Primary Text `#1D1D1F`, Secondary Text `#6E6E73`, Accent `#0071E3`, Border `#D2D2D7`
- **Typography:** `-apple-system, BlinkMacSystemFont, "SF Pro Display", "Inter", sans-serif`
- **Components:** Sticky translucent glassmorphic top navigation bar, rounded 18px cards with smooth hover elevation, `980px` pill buttons, and centered academic footer (`Created by Varun K | RA2611026050122 | Section: AIML B`) on every page.

---

## 🚀 Local Setup & Execution

```bash
git clone https://github.com/varunkolla31-create/cheminformatics-virtual-lab-apple.git
cd cheminformatics-virtual-lab-apple
pip install -r requirements.txt
streamlit run app.py
```

---

## ☁️ Deploying to Streamlit Community Cloud

1. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with GitHub.
2. Click **Create app** → **Deploy a public app from GitHub**.
3. Select repository `varunkolla31-create/cheminformatics-virtual-lab-apple`, branch `main`, and main file path `app.py`.
4. Under **Advanced settings**, choose **Python 3.11** and click **Deploy!**.

"""
Page 4: Dataset Explorer (CSV/SDF Upload, Filtering & Chemical Space Charts)
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import pandas as pd
import plotly.express as px
import streamlit as st
from rdkit import Chem
from utils import (
    apply_apple_plotly_style,
    BENCHMARK_LIBRARY,
    calculate_molecular_properties,
    init_session_state,
    inject_apple_css,
    record_studied_molecule,
    render_footer,
    render_top_nav,
)

st.set_page_config(page_title="Dataset Explorer | Cheminformatics Lab", page_icon="📈", layout="wide")
init_session_state()
inject_apple_css()
render_top_nav("Dataset Explorer")

st.markdown('<div class="apple-section-title">Chemical Dataset Explorer & SAR Analytics</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Upload a CSV or SDF dataset—or explore the built-in 24-compound pharmaceutical library—with interactive filtering and Plotly visualisations.</div>',
    unsafe_allow_html=True,
)

source_mode = st.radio(
    "Select Dataset Source",
    ["Built-in Pharmaceutical Benchmark Library (24 Compounds)", "Upload Custom CSV / SDF File", "Paste Batch SMILES Strings"],
    horizontal=True,
)

records = []

if source_mode == "Upload Custom CSV / SDF File":
    uploaded = st.file_uploader("Upload a CSV (with a 'SMILES' column) or an SDF (.sdf) file", type=["csv", "sdf"])
    if uploaded is not None:
        if uploaded.name.lower().endswith(".csv"):
            df_up = pd.read_csv(uploaded)
            smi_col = next((c for c in df_up.columns if c.lower() == "smiles"), df_up.columns[0])
            name_col = next((c for c in df_up.columns if c.lower() in ("name", "title", "id", "compound")), None)
            for idx, r in df_up.iterrows():
                cname = str(r[name_col]) if name_col else f"Compound_{idx+1}"
                p = calculate_molecular_properties(str(r[smi_col]), custom_name=cname)
                if p:
                    records.append(p)
        else:
            sdf_text = uploaded.read().decode("utf-8", errors="ignore")
            suppl = Chem.SDMolSupplier()
            suppl.SetData(sdf_text)
            for idx, mol in enumerate(suppl):
                if mol is not None:
                    smi = Chem.MolToSmiles(mol)
                    cname = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else f"SDF_Mol_{idx+1}"
                    p = calculate_molecular_properties(smi, custom_name=cname)
                    if p:
                        records.append(p)
elif source_mode == "Paste Batch SMILES Strings":
    raw = st.text_area(
        "Enter one SMILES string per line",
        value="\n".join([x["SMILES"] for x in BENCHMARK_LIBRARY[:8]]),
        height=140,
    )
    for line in raw.splitlines():
        if line.strip():
            p = calculate_molecular_properties(line.strip())
            if p:
                records.append(p)
else:
    for item in BENCHMARK_LIBRARY:
        p = calculate_molecular_properties(item["SMILES"], custom_name=item["Name"])
        if p:
            records.append(p)

if not records:
    st.info("Upload a valid CSV/SDF file or select the built-in library above to begin analysis.")
    render_footer()
    st.stop()

for r in records[:10]:
    record_studied_molecule(r)

df = pd.DataFrame(records)

st.markdown("#### Filter Chemical Space")
f1, f2, f3 = st.columns(3)
with f1:
    mw_range = st.slider("Molecular Weight Range (Da)", 0.0, float(max(600.0, df["Molecular Weight"].max())), (0.0, 560.0))
with f2:
    logp_range = st.slider("LogP Range", -5.0, 8.0, (-3.0, 6.0))
with f3:
    ro5_filter = st.selectbox("Lipinski Rule-of-Five Status", ["All Compounds", "Pass Only", "Fail Only"])

mask = (
    (df["Molecular Weight"] >= mw_range[0])
    & (df["Molecular Weight"] <= mw_range[1])
    & (df["LogP"] >= logp_range[0])
    & (df["LogP"] <= logp_range[1])
)
if ro5_filter == "Pass Only":
    mask &= df["Lipinski Check"] == "Pass"
elif ro5_filter == "Fail Only":
    mask &= df["Lipinski Check"] == "Fail"

df_filtered = df[mask].reset_index(drop=True)

c_scatter, c_bar = st.columns(2)
with c_scatter:
    fig_sc = px.scatter(
        df_filtered,
        x="Molecular Weight",
        y="LogP",
        size="TPSA",
        color="Lipinski Check",
        hover_name="Name",
        hover_data=["Formula", "SMILES", "H-bond Donors", "H-bond Acceptors"],
        color_discrete_map={"Pass": "#0071E3", "Fail": "#E03131"},
        title="Chemical Space: MW vs. LogP (Bubble Size = TPSA)",
    )
    fig_sc.add_hline(y=5.0, line_dash="dash", line_color="#E03131", annotation_text="Max LogP (5.0)")
    fig_sc.add_vline(x=500.0, line_dash="dash", line_color="#E03131", annotation_text="Max MW (500 Da)")
    apply_apple_plotly_style(fig_sc, height=400)
    st.plotly_chart(fig_sc, use_container_width=True)

with c_bar:
    metric_choice = st.selectbox(
        "Select Descriptor for Distribution Comparison",
        ["Molecular Weight", "LogP", "TPSA", "H-bond Donors", "H-bond Acceptors", "Rotatable Bonds", "Fraction Csp3"],
    )
    fig_bar = px.bar(
        df_filtered,
        x="Name",
        y=metric_choice,
        color="Lipinski Check",
        color_discrete_map={"Pass": "#0071E3", "Fail": "#E03131"},
        title=f"{metric_choice} Across Filtered Dataset",
    )
    apply_apple_plotly_style(fig_bar, height=400)
    st.plotly_chart(fig_bar, use_container_width=True)

st.markdown(f"#### Filtered Dataset Table ({len(df_filtered)} of {len(df)} compounds)")
st.dataframe(df_filtered, use_container_width=True, hide_index=True)

st.download_button(
    "Export Filtered Dataset (CSV)",
    data=df_filtered.to_csv(index=False).encode("utf-8"),
    file_name="filtered_cheminformatics_dataset.csv",
    mime="text/csv",
    type="primary",
    use_container_width=True,
)

render_footer()

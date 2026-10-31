"""
Page 2: Property Calculator (MW, LogP, TPSA, HBD/HBA, Lipinski Rule of Five)
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from utils import (
    EXAMPLE_MOLECULES,
    calculate_molecular_properties,
    generate_pdf_report,
    init_session_state,
    inject_apple_css,
    record_studied_molecule,
    render_footer,
    render_top_nav,
    validate_smiles,
)

st.set_page_config(page_title="Property Calculator | Cheminformatics Lab", page_icon="📊", layout="wide")
init_session_state()
inject_apple_css()
render_top_nav("Property Calculator")

st.markdown('<div class="apple-section-title">Physicochemical Property Calculator</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Evaluate molecular descriptors, oral bioavailability thresholds, and Lipinski\'s Rule of Five.</div>',
    unsafe_allow_html=True,
)


def _on_preset_change() -> None:
    chosen = st.session_state.get("prop_preset_select", "Custom SMILES")
    if chosen != "Custom SMILES" and chosen in EXAMPLE_MOLECULES:
        st.session_state["prop_smiles"] = EXAMPLE_MOLECULES[chosen]
        st.session_state["prop_smiles_input"] = EXAMPLE_MOLECULES[chosen]


if "prop_smiles_input" not in st.session_state:
    st.session_state["prop_smiles_input"] = st.session_state.get("prop_smiles", "CC(=O)Oc1ccccc1C(=O)O")

col_in, col_pre = st.columns([2, 1])
with col_pre:
    preset = st.selectbox(
        "Load Benchmark Molecule",
        ["Custom SMILES"] + list(EXAMPLE_MOLECULES.keys()),
        key="prop_preset_select",
        on_change=_on_preset_change,
    )

with col_in:
    smiles = st.text_input(
        "Molecule SMILES",
        key="prop_smiles_input",
        placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O",
    )
    st.session_state["prop_smiles"] = smiles

_, err = validate_smiles(smiles)
if err:
    st.error(f"Invalid SMILES: {err}")
    render_footer()
    st.stop()

with st.spinner("Calculating physicochemical properties..."):
    props = calculate_molecular_properties(smiles)

if not props:
    st.error("Unable to calculate molecular properties for the provided SMILES string.")
    render_footer()
    st.stop()

record_studied_molecule(props)

badge_cls = "apple-badge-pass" if props["Lipinski Check"] == "Pass" else "apple-badge-fail"
st.markdown(
    f"**Compound:** `{props['Name']}` (`{props['Formula']}`) &nbsp;&nbsp;|&nbsp;&nbsp; "
    f"<span class='{badge_cls}'>Lipinski Rule of Five: {props['Lipinski Check']} "
    f"({props['Lipinski Violations']} Violations)</span>",
    unsafe_allow_html=True,
)
st.write("")

r1 = st.columns(4)
r1[0].metric("Molecular Weight", f"{props['Molecular Weight']} Da", "Target <= 500 Da")
r1[1].metric("LogP (Crippen)", f"{props['LogP']}", "Target <= 5.0")
r1[2].metric("TPSA", f"{props['TPSA']} Å²", "Target < 140 Å²")
r1[3].metric("Fraction Csp3", f"{props['Fraction Csp3']}", "Saturation Index")

r2 = st.columns(4)
r2[0].metric("H-Bond Donors (HBD)", props["H-bond Donors"], "Target <= 5")
r2[1].metric("H-Bond Acceptors (HBA)", props["H-bond Acceptors"], "Target <= 10")
r2[2].metric("Rotatable Bonds", props["Rotatable Bonds"], "Target <= 10")
r2[3].metric("Ring Count", props["Ring Count"], f"{props['Heavy Atoms']} Heavy Atoms")

st.write("")
col_radar, col_table = st.columns([1, 1])

with col_radar:
    st.markdown("#### Lipinski Bioavailability Radar (Normalised to Threshold = 1.0)")
    categories = ["MW / 500", "LogP / 5", "HBD / 5", "HBA / 10", "TPSA / 140"]
    values = [
        round(props["Molecular Weight"] / 500.0, 2),
        round(max(0.0, props["LogP"]) / 5.0, 2),
        round(props["H-bond Donors"] / 5.0, 2),
        round(props["H-bond Acceptors"] / 10.0, 2),
        round(props["TPSA"] / 140.0, 2),
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Scatterpolar(
            r=[1, 1, 1, 1, 1],
            theta=categories,
            fill="toself",
            name="Ro5 Limit (1.0)",
            line_color="#D2D2D7",
        )
    )
    fig.add_trace(
        go.Scatterpolar(
            r=values,
            theta=categories,
            fill="toself",
            name=props["Name"],
            line_color="#0071E3",
        )
    )
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, max(1.5, max(values) + 0.2)])),
        margin=dict(l=30, r=30, t=30, b=30),
        height=340,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    st.plotly_chart(fig, use_container_width=True)

with col_table:
    st.markdown("#### Rule-of-Five Criterion Breakdown")
    ro5_df = pd.DataFrame([
        {
            "Rule": "Molecular Weight <= 500 Da",
            "Value": f"{props['Molecular Weight']} Da",
            "Status": "Pass" if props["Molecular Weight"] <= 500 else "Violation",
        },
        {
            "Rule": "LogP <= 5.0",
            "Value": f"{props['LogP']}",
            "Status": "Pass" if props["LogP"] <= 5.0 else "Violation",
        },
        {
            "Rule": "H-Bond Donors <= 5",
            "Value": str(props["H-bond Donors"]),
            "Status": "Pass" if props["H-bond Donors"] <= 5 else "Violation",
        },
        {
            "Rule": "H-Bond Acceptors <= 10",
            "Value": str(props["H-bond Acceptors"]),
            "Status": "Pass" if props["H-bond Acceptors"] <= 10 else "Violation",
        },
        {
            "Rule": "Veber TPSA <= 140 Å²",
            "Value": f"{props['TPSA']} Å²",
            "Status": "Pass" if props["TPSA"] <= 140 else "High Polar Area",
        },
    ])
    st.dataframe(ro5_df, use_container_width=True, hide_index=True)

    df_out = pd.DataFrame([props])
    dl1, dl2 = st.columns(2)
    with dl1:
        st.download_button(
            "Export CSV",
            data=df_out.to_csv(index=False).encode("utf-8"),
            file_name=f"properties_{props['Name']}.csv",
            mime="text/csv",
            use_container_width=True,
        )
    with dl2:
        pdf_bytes = generate_pdf_report(
            title="Molecular Properties Report",
            student_name=st.session_state["student_name"],
            reg_number=st.session_state["reg_number"],
            properties_list=[props],
        )
        st.download_button(
            "Export PDF Report",
            data=pdf_bytes,
            file_name=f"properties_{props['Name']}.pdf",
            mime="application/pdf",
            type="primary",
            use_container_width=True,
        )

render_footer()

"""
Home Page — Cheminformatics Virtual Laboratory (Apple UI/UX Edition)
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import pandas as pd
import streamlit as st
from utils import (
    BENCHMARK_LIBRARY,
    calculate_molecular_properties,
    init_session_state,
    inject_apple_css,
    render_footer,
    render_top_nav,
)

st.set_page_config(
    page_title="Cheminformatics Virtual Lab",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

init_session_state()
inject_apple_css()
render_top_nav("Home")

st.markdown(
    """
    <div class="apple-hero">
        <div class="apple-eyebrow">Computational Chemistry • Drug Discovery Suite</div>
        <div class="apple-hero-title">Cheminformatics. Designed for clarity.</div>
        <div class="apple-hero-subtitle">
            Explore 2D chemical depictions, force-field optimised 3D conformers, Lipinski drug-likeness rules,
            Tanimoto fingerprint similarity, and interactive chemical space analytics in a clean, minimalist workspace.
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

c1, c2, c3, c4, c5 = st.columns(5)
with c1:
    st.page_link("pages/1_Molecule_Viewer.py", label="🔬 Molecule Viewer", use_container_width=True)
with c2:
    st.page_link("pages/2_Property_Calculator.py", label="📊 Property Calculator", use_container_width=True)
with c3:
    st.page_link("pages/3_Similarity_Search.py", label="🧬 Similarity Search", use_container_width=True)
with c4:
    st.page_link("pages/4_Dataset_Explorer.py", label="📈 Dataset Explorer", use_container_width=True)
with c5:
    st.page_link("pages/5_About.py", label="📘 Theory & Quiz", use_container_width=True)

st.write("")
st.markdown('<div class="apple-section-title">Platform Capabilities</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Built on RDKit, py3Dmol, and Plotly with an Apple-inspired minimalist interface.</div>',
    unsafe_allow_html=True,
)

cards = [
    (
        "01 • STRUCTURE RENDERING",
        "Molecule Viewer (2D & 3D)",
        "Convert SMILES strings into crisp RDKit 2D vector depictions and interactive MMFF94 force-field optimised 3D models via py3Dmol.",
    ),
    (
        "02 • ADMET & DRUG-LIKENESS",
        "Property Calculator",
        "Compute Molecular Weight, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, Csp3 fraction, and Lipinski Rule-of-Five compliance with CSV/PDF export.",
    ),
    (
        "03 • VIRTUAL SCREENING",
        "Similarity Search",
        "Generate Morgan (ECFP4), MACCS Keys, and RDKit topological fingerprints to rank compound libraries by Tanimoto similarity.",
    ),
    (
        "04 • CHEMICAL SPACE",
        "Dataset Explorer",
        "Upload custom CSV or SDF files, filter across physicochemical thresholds, and visualise MW vs. LogP chemical space in Plotly.",
    ),
    (
        "05 • FOUNDATIONS",
        "Interactive Theory Modules",
        "Master SMILES grammar, molecular graph theory, topological polar surface area, and fingerprint mathematics through guided examples.",
    ),
    (
        "06 • EVALUATION",
        "Self-Assessment & Lab Portfolio",
        "Validate your understanding with an interactive multiple-choice quiz and export a stamped PDF laboratory portfolio.",
    ),
]

row1 = st.columns(3)
for col, (tag, title, desc) in zip(row1, cards[:3]):
    with col:
        st.markdown(
            f"""
            <div class="apple-card">
                <div class="apple-card-tag">{tag}</div>
                <h3>{title}</h3>
                <p>{desc}</p>
            </div>
            """,
            unsafe_allow_html=True,
        )

row2 = st.columns(3)
for col, (tag, title, desc) in zip(row2, cards[3:]):
    with col:
        st.markdown(
            f"""
            <div class="apple-card">
                <div class="apple-card-tag">{tag}</div>
                <h3>{title}</h3>
                <p>{desc}</p>
            </div>
            """,
            unsafe_allow_html=True,
        )

st.write("")
st.markdown('<div class="apple-section-title">Benchmark Reference Compounds</div>', unsafe_allow_html=True)
st.caption("Use any of these curated pharmaceutical and organic structures across all laboratory modules.")

bench_rows = []
for item in BENCHMARK_LIBRARY[:10]:
    p = calculate_molecular_properties(item["SMILES"], custom_name=item["Name"])
    if p:
        bench_rows.append(p)

df_bench = pd.DataFrame(bench_rows)[
    ["Name", "Formula", "SMILES", "Molecular Weight", "LogP", "TPSA", "Lipinski Check"]
]
st.dataframe(df_bench, use_container_width=True, hide_index=True)

render_footer()

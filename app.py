"""
Home Page — Cheminformatics Virtual Laboratory
SRM Institute of Science and Technology
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import pandas as pd
import streamlit as st
from utils import (
    BENCHMARK_LIBRARY,
    calculate_molecular_properties,
    init_session_state,
    inject_apple_css,
    render_2d_image,
    render_footer,
    render_top_nav,
)

st.set_page_config(
    page_title="Cheminformatics Virtual Lab | SRM IST",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

init_session_state()
inject_apple_css()
render_top_nav("Home")

# Hero Section — Clean, Asymmetric Product Header
st.markdown(
    """
    <div class="apple-hero">
        <div class="apple-eyebrow">
            SRM Institute of Science and Technology, Tiruchirappalli
        </div>
        <div style="font-size: 0.86rem; font-weight: 500; color: #6E6E73; margin-bottom: 0.9rem;">
            Faculty of Engineering and Technology &nbsp;•&nbsp; School of Sciences &nbsp;•&nbsp; Division of Chemistry
        </div>
        <div class="apple-hero-title">Cheminformatics Virtual Laboratory</div>
        <div class="apple-hero-subtitle">
            An interactive computational chemistry environment for molecular structure representation,
            3D conformer geometry optimisation, physicochemical descriptor analysis, and virtual screening.
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

# Direct Module Navigation Bar
c1, c2, c3, c4, c5 = st.columns(5)
with c1:
    st.page_link("pages/1_Molecule_Viewer.py", label="Molecule Viewer", use_container_width=True)
with c2:
    st.page_link("pages/2_Property_Calculator.py", label="Property Calculator", use_container_width=True)
with c3:
    st.page_link("pages/3_Similarity_Search.py", label="Similarity Search", use_container_width=True)
with c4:
    st.page_link("pages/4_Dataset_Explorer.py", label="Dataset Explorer", use_container_width=True)
with c5:
    st.page_link("pages/5_About.py", label="Theory & Portfolio", use_container_width=True)

st.write("")

# Core Laboratory Modules Grid
st.markdown('<div class="apple-section-title">Laboratory Workbench</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Select a tool below or use the sidebar to move through the experimental workflow.</div>',
    unsafe_allow_html=True,
)

cards = [
    (
        "Structure & Conformation",
        "2D & 3D Molecule Viewer",
        "Parse SMILES strings into canonical RDKit 2D depictions and interactive MMFF94 force-field minimised 3D conformers with SDF and PNG export.",
    ),
    (
        "Physicochemical Profiling",
        "Molecular Property Calculator",
        "Calculate Molecular Weight, Crippen LogP, TPSA, hydrogen-bond counts, rotatable bonds, and Lipinski Rule-of-Five compliance.",
    ),
    (
        "Lead Identification",
        "Fingerprint Similarity Screening",
        "Generate Morgan (ECFP4), MACCS Keys, and topological bit vectors to rank reference compounds by Tanimoto similarity.",
    ),
    (
        "Structure–Property Analysis",
        "Dataset Explorer & Chemical Space",
        "Upload custom CSV or SDF compound libraries, apply multi-parameter filters, and inspect MW vs. LogP distributions.",
    ),
    (
        "Foundational Concepts",
        "Theory & Self-Assessment",
        "Review concise notes on SMILES grammar, molecular graph representations, and oral bioavailability thresholds with a 10-question quiz.",
    ),
    (
        "Lab Documentation",
        "Session Portfolio & PDF Report",
        "Compile all molecules analysed during your session alongside your quiz score into a formatted laboratory PDF report.",
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

# Reference Molecules Gallery & Benchmark Table
st.markdown('<div class="apple-section-title">Reference Structures</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Standard organic and pharmaceutical compounds pre-loaded for immediate testing.</div>',
    unsafe_allow_html=True,
)

featured_names = ["Aspirin", "Caffeine", "Paracetamol", "Ibuprofen", "Benzene", "Ethanol"]
featured_items = [x for x in BENCHMARK_LIBRARY if x["Name"] in featured_names]

mol_cols = st.columns(6)
for col, item in zip(mol_cols, featured_items):
    with col:
        img = render_2d_image(item["SMILES"], width=240, height=180)
        if img:
            st.image(img, use_container_width=True)
        st.markdown(
            f"<div style='text-align:center; font-size:0.84rem;'><strong>{item['Name']}</strong><br/>"
            f"<code style='font-size:0.72rem; color:#6E6E73;'>{item['SMILES']}</code></div>",
            unsafe_allow_html=True,
        )

st.write("")
bench_rows = []
for item in BENCHMARK_LIBRARY[:10]:
    p = calculate_molecular_properties(item["SMILES"], custom_name=item["Name"])
    if p:
        bench_rows.append(p)

df_bench = pd.DataFrame(bench_rows)[
    ["Name", "Formula", "SMILES", "Molecular Weight", "LogP", "TPSA", "Lipinski Check"]
]
st.dataframe(df_bench, use_container_width=True, hide_index=True)

st.write("")

# Clean Academic Attribution Cards
fac_col, stu_col = st.columns(2)

with fac_col:
    st.markdown(
        """
        <div class="apple-card">
            <div class="apple-card-tag">Course Coordinator &amp; Faculty Mentor</div>
            <h3>Dr. P Pachamuthu</h3>
            <p style="margin-bottom:0.55rem;">
                Professor, Division of Chemistry, School of Sciences<br/>
                Faculty of Engineering and Technology<br/>
                SRM Institute of Science and Technology, Tiruchirappalli
            </p>
            <p style="font-size:0.85rem;">
                Contact: <a href="mailto:pachamuthup@srmist.edu.in" style="color:#0071E3; text-decoration:none;">pachamuthup@srmist.edu.in</a>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

with stu_col:
    st.markdown(
        """
        <div class="apple-card">
            <div class="apple-card-tag">Student in SRMIST</div>
            <h3>Varun K</h3>
            <p style="margin-bottom:0.55rem;">
                Registration No: <strong>RA2611026050122</strong> &nbsp;•&nbsp; Section: <strong>AIML B</strong><br/>
                B.Tech in Artificial Intelligence and Machine Learning<br/>
                SRM Institute of Science and Technology
            </p>
            <p style="font-size:0.85rem; color:#6E6E73;">
                Built with Python 3.11, Streamlit, RDKit, py3Dmol &amp; Plotly
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

render_footer()

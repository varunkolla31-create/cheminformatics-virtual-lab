"""
Home Page — Cheminformatics Virtual Laboratory
Modelled on SRM Institute of Science and Technology Cheminformatics Virtual Laboratory
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
        <div class="apple-eyebrow">
            SRM Institute of Science and Technology, Tiruchirappalli
        </div>
        <div style="font-size: 0.92rem; font-weight: 500; color: #6E6E73; margin-bottom: 0.85rem;">
            Faculty of Engineering and Technology &nbsp;|&nbsp; School of Sciences &nbsp;|&nbsp; Division of Chemistry
        </div>
        <div class="apple-hero-title">🧪 Cheminformatics Virtual Laboratory</div>
        <div class="apple-hero-subtitle">
            This interactive virtual laboratory helps students understand important concepts in
            cheminformatics and computational chemistry—reimagined with Apple-inspired clarity.
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

c1, c2, c3, c4, c5 = st.columns(5)
with c1:
    st.page_link("pages/1_Molecule_Viewer.py", label="🧬 Molecular Visualization", use_container_width=True)
with c2:
    st.page_link("pages/2_Property_Calculator.py", label="📊 Properties Report", use_container_width=True)
with c3:
    st.page_link("pages/3_Similarity_Search.py", label="🔍 Similarity Search", use_container_width=True)
with c4:
    st.page_link("pages/4_Dataset_Explorer.py", label="📈 Structure–Property", use_container_width=True)
with c5:
    st.page_link("pages/5_About.py", label="📝 Theory & Assessment", use_container_width=True)

st.write("")

col_welcome, col_checklist = st.columns([1.15, 1.0])

with col_welcome:
    st.markdown(
        """
        <div class="apple-card">
            <div class="apple-card-tag">WELCOME TO THE VIRTUAL LABORATORY</div>
            <h3>Interactive Computational Chemistry Workspace</h3>
            <p style="margin-bottom: 1rem;">
                Designed for undergraduate engineering and science students at SRM Institute of Science and Technology,
                this platform bridges chemical structure theory with hands-on computational analysis using Python,
                Streamlit, RDKit, py3Dmol, and Plotly.
            </p>
            <div style="border-top: 1px solid #D2D2D7; padding-top: 0.9rem; margin-top: 0.9rem;">
                <strong style="color:#1D1D1F; font-size:0.9rem;">🎓 Core Learning Objectives:</strong>
                <ul style="color:#6E6E73; font-size:0.9rem; line-height:1.65; margin-top:0.4rem; margin-bottom:0;">
                    <li>SMILES notation &amp; canonical molecular representation</li>
                    <li>2D chemical structure depiction &amp; validation</li>
                    <li>Interactive 3D conformer embedding &amp; force-field optimisation</li>
                    <li>Molecular descriptors &amp; Lipinski's Rule of Five</li>
                    <li>Structure–property relationships (SAR) &amp; Tanimoto similarity</li>
                </ul>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col_checklist:
    st.markdown(
        """
        <div class="apple-card">
            <div class="apple-card-tag">STUDENT CAPABILITIES</div>
            <h3>🧪 What can students do?</h3>
            <div style="display: grid; grid-template-columns: 1fr; gap: 0.48rem; margin-top: 0.75rem; font-size: 0.93rem; color: #1D1D1F;">
                <div>✅ <strong>Enter any SMILES notation</strong> and validate molecular structures</div>
                <div>✅ <strong>Generate 2D structures</strong> with publication-quality RDKit rendering</div>
                <div>✅ <strong>Generate interactive 3D structures</strong> with MMFF94 geometry optimisation</div>
                <div>✅ <strong>Select different molecular visualization styles</strong> (Stick, Sphere, Line, CPK)</div>
                <div>✅ <strong>Calculate complete molecular properties</strong> (MW, LogP, TPSA, HBD, HBA, Rings)</div>
                <div>✅ <strong>Download molecular reports</strong> in formatted CSV and official PDF formats</div>
                <div>✅ <strong>Compare multiple molecules</strong> across chemical space &amp; Tanimoto similarity</div>
                <div>✅ <strong>Complete an interactive assessment</strong> and export a final project report</div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.write("")

st.markdown('<div class="apple-section-title">Virtual Laboratory Modules</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Navigate through the guided workflow from foundational theory to your final downloadable project report.</div>',
    unsafe_allow_html=True,
)

cards = [
    (
        "01 • MOLECULAR VISUALIZATION",
        "🧬 2D & 3D Structure Viewer",
        "Enter any SMILES string, validate chemical syntax, inspect 2D depictions, and rotate interactive 3D conformers in stick, sphere, or line styles.",
    ),
    (
        "02 • MOLECULAR PROPERTIES REPORT",
        "📊 Physicochemical Descriptors",
        "Compute Molecular Weight, Crippen LogP, TPSA, H-bond donors/acceptors, rotatable bonds, and Lipinski Rule-of-Five status with instant CSV & PDF export.",
    ),
    (
        "03 • SIMILARITY SEARCH",
        "🔍 Fingerprint Virtual Screening",
        "Compute Morgan (ECFP4), MACCS Keys, and RDKit topological fingerprints to rank reference molecules by Tanimoto similarity.",
    ),
    (
        "04 • STRUCTURE–PROPERTY ANALYSIS",
        "📈 Multi-Molecule Comparison",
        "Compare multiple molecules side-by-side or upload CSV/SDF datasets to explore MW vs. LogP chemical space and descriptor distributions.",
    ),
    (
        "05 • THEORY & ASSESSMENT",
        "📝 Self-Assessment Quiz",
        "Study concise modules on SMILES, molecular graphs, and drug-likeness rules, then complete the 10-question graded assessment.",
    ),
    (
        "06 • FINAL PROJECT REPORT",
        "📄 Downloadable Lab Portfolio",
        "Compile your student registration profile, studied molecules table, summary statistics, and assessment score into a single academic PDF.",
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

st.markdown('<div class="apple-section-title">Example Molecules</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Explore these benchmark organic and pharmaceutical structures across the Virtual Laboratory.</div>',
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
            f"<div style='text-align:center; font-size:0.88rem;'><strong>{item['Name']}</strong><br/>"
            f"<code style='font-size:0.75rem;'>{item['SMILES']}</code></div>",
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

fac_col, stu_col = st.columns(2)

with fac_col:
    st.markdown(
        """
        <div class="apple-card">
            <div class="apple-card-tag">FACULTY REFERENCE &amp; GUIDANCE</div>
            <h3>👨‍🏫 Developed by</h3>
            <p style="color:#1D1D1F; font-weight:600; font-size:1.02rem; margin-bottom:0.25rem;">
                Dr. P Pachamuthu
            </p>
            <p style="margin-bottom:0.6rem;">
                Professor<br/>
                Division of Chemistry, School of Sciences<br/>
                Faculty of Engineering and Technology<br/>
                SRM Institute of Science and Technology, Tiruchirappalli
            </p>
            <p style="font-size:0.88rem;">
                📧 <strong>For queries, contact:</strong>
                <a href="mailto:pachamuthup@srmist.edu.in" style="color:#0071E3; text-decoration:none;">pachamuthup@srmist.edu.in</a>
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

with stu_col:
    st.markdown(
        """
        <div class="apple-card">
            <div class="apple-card-tag">STUDENT IN SRMIST</div>
            <h3>🎓 Student Developer Profile</h3>
            <p style="color:#1D1D1F; font-weight:600; font-size:1.02rem; margin-bottom:0.25rem;">
                Varun K
            </p>
            <p style="margin-bottom:0.6rem;">
                Registration Number: <strong>RA2611026050122</strong><br/>
                Section: <strong>AIML B</strong> (B.Tech — Artificial Intelligence &amp; Machine Learning)<br/>
                SRM Institute of Science and Technology
            </p>
            <p style="font-size:0.88rem; color:#6E6E73;">
                🧪 Cheminformatics Virtual Laboratory &nbsp;|&nbsp; Developed using Python 3.11, Streamlit, RDKit, py3Dmol &amp; Plotly
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

render_footer()

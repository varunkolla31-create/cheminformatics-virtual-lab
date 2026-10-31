"""
Page 5: About, Theory Foundations, Self-Assessment Quiz & Final Portfolio
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import pandas as pd
import streamlit as st
from utils import (
    BENCHMARK_LIBRARY,
    QUIZ_QUESTIONS,
    calculate_molecular_properties,
    generate_pdf_report,
    init_session_state,
    inject_apple_css,
    render_footer,
    render_top_nav,
)

st.set_page_config(page_title="About & Theory | Cheminformatics Lab", page_icon="📘", layout="wide")
init_session_state()
inject_apple_css()
render_top_nav("About")

st.markdown('<div class="apple-section-title">Cheminformatics Theory, Assessment & Portfolio</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Review foundational concepts, complete the 10-question self-assessment quiz, and download your stamped PDF lab portfolio.</div>',
    unsafe_allow_html=True,
)

tab_theory, tab_quiz, tab_portfolio = st.tabs([
    "Core Theory",
    "Self-Assessment Quiz",
    "Final Lab Portfolio",
])

with tab_theory:
    with st.expander("1. SMILES Notation & Canonicalisation", expanded=True):
        st.markdown(
            """
            - **SMILES (Simplified Molecular Input Line Entry System)** encodes molecular graphs as compact ASCII text.
            - **Atoms & Bonds:** Heavy atoms use atomic symbols (`C`, `N`, `O`, `F`, `Cl`); single bonds are implicit, double bonds use `=`, triple bonds use `#`.
            - **Branches & Rings:** Parentheses denote side chains (`CC(=O)O` for acetic acid) and matching digits close rings (`c1ccccc1` for aromatic benzene).
            """
        )
    with st.expander("2. Molecular Graphs & 3D Force Fields (MMFF94)"):
        st.markdown(
            """
            - A molecule is stored computationally as a graph $G = (V, E)$ where vertices $V$ are atoms and edges $E$ are bonds.
            - 3D distance-geometry embedding (`ETKDGv3`) followed by **MMFF94** (Merck Molecular Force Field) energy minimisation yields realistic bond lengths, angles, and dihedral conformations.
            """
        )
    with st.expander("3. Lipinski's Rule of Five & TPSA"):
        st.markdown(
            """
            - **Lipinski's Rule of Five (Ro5)** predicts oral bioavailability: $\\text{MW} \\le 500\\text{ Da}$, $\\text{LogP} \\le 5$, $\\text{HBD} \\le 5$, and $\\text{HBA} \\le 10$ (at most 1 violation allowed).
            - **TPSA (Topological Polar Surface Area)** sums polar surface contributions of N and O atoms; $\\text{TPSA} < 140\\text{ \\AA}^2$ correlates with high intestinal absorption.
            """
        )
    with st.expander("4. Molecular Fingerprints & Tanimoto Similarity"):
        st.markdown(
            """
            - Circular fingerprints (**Morgan / ECFP4**) encode atom neighbourhoods into a fixed-length bit vector (e.g. 2048 bits).
            - **Tanimoto Similarity:** $T(A, B) = \\frac{|A \\cap B|}{|A| + |B| - |A \\cap B|}$ quantifies structural similarity from `0.0` (disjoint) to `1.0` (identical).
            """
        )

with tab_quiz:
    st.markdown("#### Interactive Cheminformatics Quiz")
    score = 0
    for idx, q in enumerate(QUIZ_QUESTIONS):
        st.markdown(f"**Q{idx+1}. {q['q']}**")
        choice = st.radio(
            f"Select answer for Q{idx+1}",
            q["options"],
            key=f"quiz_q_{idx}",
            label_visibility="collapsed",
        )
        st.session_state["quiz_answers"][idx] = choice
        if choice == q["answer"]:
            score += 1
        st.write("")

    if st.button("Submit & Grade Assessment", type="primary"):
        st.success(f"🏆 **Your Quiz Score: {score} / {len(QUIZ_QUESTIONS)} ({score / len(QUIZ_QUESTIONS) * 100:.0f}%)**")
        for idx, q in enumerate(QUIZ_QUESTIONS):
            user_ans = st.session_state["quiz_answers"].get(idx)
            if user_ans == q["answer"]:
                st.markdown(f"- **Q{idx+1}: Correct!** {q['explanation']}")
            else:
                st.markdown(f"- **Q{idx+1}: Incorrect** (Correct: *{q['answer']}*). {q['explanation']}")

with tab_portfolio:
    st.markdown("#### Session Portfolio & Official PDF Export")
    c1, c2, c3 = st.columns(3)
    c1.metric("Student Name", st.session_state["student_name"])
    c2.metric("Registration Number", st.session_state["reg_number"])
    c3.metric("Section", st.session_state["section"])

    studied = list(st.session_state["studied_molecules"].values())
    if not studied:
        studied = [calculate_molecular_properties(x["SMILES"], custom_name=x["Name"]) for x in BENCHMARK_LIBRARY[:5]]

    df_port = pd.DataFrame(studied)
    st.dataframe(df_port, use_container_width=True, hide_index=True)

    pdf_bytes = generate_pdf_report(
        title=f"Final Lab Portfolio ({st.session_state['section']})",
        student_name=st.session_state["student_name"],
        reg_number=st.session_state["reg_number"],
        properties_list=studied,
        assessment_score=(score, len(QUIZ_QUESTIONS)),
    )
    st.download_button(
        "Download Lab Portfolio (PDF)",
        data=pdf_bytes,
        file_name=f"Cheminformatics_Portfolio_{st.session_state['reg_number']}.pdf",
        mime="application/pdf",
        type="primary",
        use_container_width=True,
    )

render_footer()

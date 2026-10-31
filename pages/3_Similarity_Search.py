"""
Page 3: Similarity Search (Tanimoto Similarity & Molecular Fingerprints)
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import plotly.express as px
import streamlit as st
from utils import (
    EXAMPLE_MOLECULES,
    compute_similarity_search,
    init_session_state,
    inject_apple_css,
    render_2d_image,
    render_footer,
    render_top_nav,
    validate_smiles,
)

st.set_page_config(page_title="Similarity Search | Cheminformatics Lab", page_icon="🧬", layout="wide")
init_session_state()
inject_apple_css()
render_top_nav("Similarity Search")

st.markdown('<div class="apple-section-title">🧬 Molecular Fingerprint & Tanimoto Similarity Search</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Screen a curated reference library against your query compound using 2D molecular fingerprints.</div>',
    unsafe_allow_html=True,
)

c1, c2, c3, c4 = st.columns([2.2, 1.2, 1.2, 1.0])
with c1:
    query_smiles = st.text_input("Query Molecule SMILES", value=st.session_state["sim_smiles"])
    st.session_state["sim_smiles"] = query_smiles
with c2:
    preset = st.selectbox("Quick Query Preset", ["Custom"] + list(EXAMPLE_MOLECULES.keys()))
    if preset != "Custom":
        st.session_state["sim_smiles"] = EXAMPLE_MOLECULES[preset]
        query_smiles = st.session_state["sim_smiles"]
with c3:
    fp_type = st.selectbox("Fingerprint Algorithm", ["Morgan (ECFP4)", "MACCS Keys (166-bit)", "RDKit Topological"])
with c4:
    min_sim = st.slider("Min Tanimoto Cutoff", 0.0, 1.0, 0.10, 0.05)

_, err = validate_smiles(query_smiles)
if err:
    st.error(f"⚠️ {err}")
    render_footer()
    st.stop()

with st.spinner("Computing molecular fingerprints and Tanimoto coefficients..."):
    df_sim = compute_similarity_search(query_smiles, fp_type=fp_type)

df_filtered = df_sim[df_sim["Tanimoto Similarity"] >= min_sim].reset_index(drop=True)

top_cols = st.columns(4)
for col, (_, row) in zip(top_cols, df_sim.head(4).iterrows()):
    with col:
        img = render_2d_image(row["SMILES"], width=300, height=220)
        if img:
            st.image(img, use_container_width=True)
        st.markdown(
            f"**{row['Name']}**  \n`Tanimoto = {row['Tanimoto Similarity']:.3f}` • {row['Category']}"
        )

st.divider()
left_chart, right_table = st.columns([1.1, 1.2])

with left_chart:
    fig = px.bar(
        df_filtered.head(12),
        x="Tanimoto Similarity",
        y="Name",
        orientation="h",
        color="Tanimoto Similarity",
        color_continuous_scale=["#D2D2D7", "#0071E3"],
        title=f"Top Ranked Analogs ({fp_type})",
    )
    fig.update_layout(yaxis={"categoryorder": "total ascending"}, height=420)
    st.plotly_chart(fig, use_container_width=True)

with right_table:
    st.markdown(f"#### Screening Hits ({len(df_filtered)} compounds >= {min_sim:.2f})")
    st.dataframe(df_filtered, use_container_width=True, hide_index=True, height=360)
    st.download_button(
        "📥 Download Similarity Ranking (CSV)",
        data=df_filtered.to_csv(index=False).encode("utf-8"),
        file_name="tanimoto_similarity_results.csv",
        mime="text/csv",
        type="primary",
        use_container_width=True,
    )

render_footer()

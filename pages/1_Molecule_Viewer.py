"""
Page 1: Molecule Viewer (2D RDKit + 3D py3Dmol)
Created by Varun K | RA2611026050122 | Section: AIML B
"""

import py3Dmol
import streamlit as st
from utils import (
    EXAMPLE_MOLECULES,
    build_3d_molblock,
    calculate_molecular_properties,
    init_session_state,
    inject_apple_css,
    record_studied_molecule,
    render_2d_image,
    render_3d_viewer,
    render_footer,
    render_top_nav,
    validate_smiles,
)

st.set_page_config(page_title="Molecule Viewer | Cheminformatics Lab", page_icon="🔬", layout="wide")
init_session_state()
inject_apple_css()
render_top_nav("Molecule Viewer")

st.markdown('<div class="apple-section-title">🔬 Molecule Viewer (2D & 3D)</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="apple-section-lead">Inspect publication-grade 2D skeletal structures and MMFF94 force-field minimised 3D conformers.</div>',
    unsafe_allow_html=True,
)

st.caption("Quick-Load Benchmark Molecules")
cols = st.columns(len(EXAMPLE_MOLECULES))
for col, (name, smi) in zip(cols, EXAMPLE_MOLECULES.items()):
    with col:
        if st.button(name, key=f"btn_{name}", use_container_width=True):
            st.session_state["viz_smiles"] = smi

smiles_input = st.text_input(
    "SMILES Notation Input",
    value=st.session_state["viz_smiles"],
    placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O",
)
st.session_state["viz_smiles"] = smiles_input

mol, err = validate_smiles(smiles_input)
if err:
    st.error(f"⚠️ {err}")
    render_footer()
    st.stop()

props = calculate_molecular_properties(smiles_input)
if props:
    record_studied_molecule(props)

view_mode = st.radio(
    "Select Visualisation Mode",
    options=["Side-by-Side (2D & 3D)", "2D Skeletal Depiction Only", "3D Interactive Conformer Only"],
    horizontal=True,
)

col_style, col_scheme, col_bg = st.columns(3)
with col_style:
    style_3d = st.selectbox("3D Style", ["stick", "sphere", "line", "cross"], index=0)
with col_scheme:
    color_scheme = st.selectbox(
        "Atom Colour Scheme",
        ["default (CPK)", "cyanCarbon", "greenCarbon", "magentaCarbon", "yellowCarbon"],
        index=0,
    )
with col_bg:
    bg_choice = st.selectbox("3D Canvas Background", ["Pure White (#FFFFFF)", "Apple Light (#F5F5F7)", "Dark Slate (#1D1D1F)"])

bg_hex = "#FFFFFF" if "White" in bg_choice else ("#F5F5F7" if "Light" in bg_choice else "#1D1D1F")

img_bytes = None
molblock = None

if view_mode in ("Side-by-Side (2D & 3D)", "2D Skeletal Depiction Only"):
    with st.spinner("Rendering 2D chemical structure..."):
        img_bytes = render_2d_image(smiles_input, width=540, height=400)

if view_mode in ("Side-by-Side (2D & 3D)", "3D Interactive Conformer Only"):
    with st.spinner("Optimising 3D geometry with MMFF94 force field..."):
        molblock, err_3d = build_3d_molblock(smiles_input)
        if err_3d:
            st.warning(err_3d)

if view_mode == "Side-by-Side (2D & 3D)":
    left, right = st.columns(2)
    with left:
        st.markdown("#### 2D Skeletal Structure (RDKit)")
        if img_bytes:
            st.image(img_bytes, use_container_width=True)
            st.download_button(
                "📥 Download 2D PNG",
                data=img_bytes,
                file_name=f"{props['Name']}_2D.png",
                mime="image/png",
                use_container_width=True,
            )
    with right:
        st.markdown("#### 3D Minimised Conformer (py3Dmol)")
        if molblock:
            scheme = "default" if color_scheme.startswith("default") else color_scheme
            view = py3Dmol.view(width=540, height=390)
            view.addModel(molblock, "mol")
            view.setStyle({style_3d: {"colorscheme": scheme}})
            view.setBackgroundColor(bg_hex)
            view.zoomTo()
            render_3d_viewer(view, height=400, width=540)
            st.download_button(
                "📥 Download 3D SDF / MolBlock",
                data=molblock.encode("utf-8"),
                file_name=f"{props['Name']}_3D.sdf",
                mime="chemical/x-mdl-sdfile",
                use_container_width=True,
            )
elif view_mode == "2D Skeletal Depiction Only":
    if img_bytes:
        st.image(img_bytes, width=560)
else:
    if molblock:
        scheme = "default" if color_scheme.startswith("default") else color_scheme
        view = py3Dmol.view(width=760, height=460)
        view.addModel(molblock, "mol")
        view.setStyle({style_3d: {"colorscheme": scheme}})
        view.setBackgroundColor(bg_hex)
        view.zoomTo()
        render_3d_viewer(view, height=470, width=760)

st.divider()
m1, m2, m3, m4 = st.columns(4)
m1.metric("Compound Name", props["Name"])
m2.metric("Molecular Formula", props["Formula"])
m3.metric("Molecular Weight", f"{props['Molecular Weight']} g/mol")
m4.metric("Heavy Atom Count", props["Heavy Atoms"])

render_footer()

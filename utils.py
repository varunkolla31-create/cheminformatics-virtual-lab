"""
Shared Cheminformatics Engine & Apple Design System Utilities
Created by Varun K | RA2611026050122 | Section: AIML B
"""

from __future__ import annotations

import io
import math
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

try:
    from rdkit.Chem import Draw as RDKitDraw
    HAS_RDKIT_DRAW = True
except Exception:
    RDKitDraw = None
    HAS_RDKIT_DRAW = False

try:
    from stmol import showmol as _stmol_showmol
    HAS_STMOL = True
except Exception:
    _stmol_showmol = None
    HAS_STMOL = False

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

# -----------------------------------------------------------------------------
# CONSTANTS & BENCHMARK MOLECULE LIBRARY
# -----------------------------------------------------------------------------

FOOTER_TEXT = "Created by Varun K | RA2611026050122 | Section: AIML B"

EXAMPLE_MOLECULES: Dict[str, str] = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Phenol": "Oc1ccccc1",
    "Acetic Acid": "CC(=O)O",
}

BENCHMARK_LIBRARY: List[Dict[str, str]] = [
    {"Name": "Aspirin", "Category": "NSAID Analgesic", "SMILES": "CC(=O)Oc1ccccc1C(=O)O"},
    {"Name": "Salicylic Acid", "Category": "Metabolite / Keratolytic", "SMILES": "O=C(O)c1ccccc1O"},
    {"Name": "Methyl Salicylate", "Category": "Topical Analgesic", "SMILES": "COC(=O)c1ccccc1O"},
    {"Name": "Paracetamol", "Category": "Analgesic / Antipyretic", "SMILES": "CC(=O)Nc1ccc(O)cc1"},
    {"Name": "Phenacetin", "Category": "Arylacetamide Analgesic", "SMILES": "CCOc1ccc(NC(C)=O)cc1"},
    {"Name": "Ibuprofen", "Category": "Propionic Acid NSAID", "SMILES": "CC(C)Cc1ccc(C(C)C(=O)O)cc1"},
    {"Name": "Naproxen", "Category": "Propionic Acid NSAID", "SMILES": "COc1ccc2cc(C(C)C(=O)O)ccc2c1"},
    {"Name": "Ketoprofen", "Category": "Propionic Acid NSAID", "SMILES": "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1"},
    {"Name": "Caffeine", "Category": "Methylxanthine Stimulant", "SMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
    {"Name": "Theobromine", "Category": "Methylxanthine Alkaloid", "SMILES": "Cn1cnc2c1c(=O)[nH]c(=O)n2C"},
    {"Name": "Theophylline", "Category": "Bronchodilator Xanthine", "SMILES": "Cn1c(=O)c2[nH]cnc2n(C)c1=O"},
    {"Name": "Dopamine", "Category": "Catecholamine Neurotransmitter", "SMILES": "NCCc1ccc(O)c(O)c1"},
    {"Name": "Adrenaline", "Category": "Hormone / Neurotransmitter", "SMILES": "CNC[C@H](O)c1ccc(O)c(O)c1"},
    {"Name": "Serotonin", "Category": "Indoleamine Neurotransmitter", "SMILES": "NCCc1c[nH]c2ccc(O)cc12"},
    {"Name": "Penicillin G", "Category": "Beta-Lactam Antibiotic", "SMILES": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"},
    {"Name": "Ciprofloxacin", "Category": "Fluoroquinolone Antibiotic", "SMILES": "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O"},
    {"Name": "Metformin", "Category": "Biguanide Antidiabetic", "SMILES": "CN(C)C(=N)NC(=N)N"},
    {"Name": "Atorvastatin", "Category": "Statin Lipid-Lowering", "SMILES": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"},
    {"Name": "Omeprazole", "Category": "Proton Pump Inhibitor", "SMILES": "COc1ccc2nc(S(=O)Cc3ncc(C)c(OC)c3C)[nH]c2c1"},
    {"Name": "Loratadine", "Category": "Tricyclic Antihistamine", "SMILES": "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1"},
    {"Name": "Ethanol", "Category": "Aliphatic Alcohol", "SMILES": "CCO"},
    {"Name": "Benzene", "Category": "Aromatic Hydrocarbon", "SMILES": "c1ccccc1"},
    {"Name": "Phenol", "Category": "Aromatic Alcohol", "SMILES": "Oc1ccccc1"},
    {"Name": "Acetic Acid", "Category": "Carboxylic Acid", "SMILES": "CC(=O)O"},
]

LEARNING_OBJECTIVES: List[str] = [
    "SMILES notation & canonicalisation",
    "Molecular graph representation",
    "2D chemical structure depiction",
    "3D conformer embedding & MMFF94 force-field optimisation",
    "Physicochemical descriptors & Lipinski's Rule of Five",
    "Molecular fingerprints & Tanimoto similarity search",
    "Structure-property relationships & dataset exploration",
]

QUIZ_QUESTIONS: List[Dict[str, Any]] = [
    {
        "q": "What does the acronym SMILES stand for in cheminformatics?",
        "options": [
            "Simplified Molecular Input Line Entry System",
            "Standard Molecular Index for Linear Element Strings",
            "Structural Modeling of Inorganic and Linear Element Systems",
            "Systematic Molecular Identification and Labeling Engine Syntax",
        ],
        "answer": "Simplified Molecular Input Line Entry System",
        "explanation": "SMILES (Simplified Molecular Input Line Entry System) represents chemical structures using short ASCII strings.",
    },
    {
        "q": "Which SMILES string correctly represents Ethanol?",
        "options": ["CCO", "CC(=O)O", "c1ccccc1", "COC"],
        "answer": "CCO",
        "explanation": "Ethanol (CH3-CH2-OH) has two aliphatic carbons single-bonded to an oxygen atom, written as CCO (hydrogens are implicit).",
    },
    {
        "q": "In SMILES notation, what do lowercase letters such as 'c' and 'n' denote?",
        "options": [
            "Atoms in an aromatic ring",
            "Chiral centers",
            "Isotopically labeled atoms",
            "Charged radical species",
        ],
        "answer": "Atoms in an aromatic ring",
        "explanation": "Lowercase atomic symbols represent aromatic atoms (e.g., 'c1ccccc1' for aromatic benzene vs. 'C1CCCCC1' for cyclohexane).",
    },
    {
        "q": "In a molecular graph representation, what do the nodes (vertices) and edges represent?",
        "options": [
            "Nodes represent atoms and edges represent chemical bonds",
            "Nodes represent bonds and edges represent lone pairs",
            "Nodes represent functional groups and edges represent dipole moments",
            "Nodes represent orbitals and edges represent hybridization states",
        ],
        "answer": "Nodes represent atoms and edges represent chemical bonds",
        "explanation": "A molecular graph models a molecule mathematically with atoms as vertices (nodes) and covalent bonds as edges.",
    },
    {
        "q": "What does the descriptor LogP (octanol-water partition coefficient) measure?",
        "options": [
            "Lipophilicity (hydrophobicity) of a molecule",
            "Molecular polarizability in vacuum",
            "Total number of rotatable bonds",
            "Acid dissociation constant in aqueous buffer",
        ],
        "answer": "Lipophilicity (hydrophobicity) of a molecule",
        "explanation": "LogP is the base-10 logarithm of a compound's concentration ratio between 1-octanol (lipid-like) and water, measuring lipophilicity.",
    },
    {
        "q": "According to Lipinski's Rule of Five, an orally active drug candidate should generally have a Molecular Weight (MW) of:",
        "options": ["<= 500 Da", "<= 100 Da", ">= 500 Da", "Exactly 250 Da"],
        "answer": "<= 500 Da",
        "explanation": "Lipinski's criteria stipulate MW <= 500 Da, LogP <= 5, H-bond donors <= 5, and H-bond acceptors <= 10 for good oral bioavailability.",
    },
    {
        "q": "What does a Tanimoto similarity coefficient of 1.0 between two molecular bit fingerprints indicate?",
        "options": [
            "Identical fingerprint bit patterns under the chosen representation",
            "Completely dissimilar chemical structures",
            "Opposite stereochemical configurations only",
            "A molecular weight difference of exactly 100 Da",
        ],
        "answer": "Identical fingerprint bit patterns under the chosen representation",
        "explanation": "Tanimoto coefficient T = c / (a + b - c) ranges from 0.0 (no shared bits) to 1.0 (identical active bits).",
    },
    {
        "q": "What property does TPSA (Topological Polar Surface Area) help predict?",
        "options": [
            "Cell membrane permeability and blood-brain barrier penetration",
            "Covalent bond dissociation energy",
            "Nuclear magnetic resonance chemical shifts",
            "Optical rotation angle of enantiomers",
        ],
        "answer": "Cell membrane permeability and blood-brain barrier penetration",
        "explanation": "TPSA sums the surface contributions of polar atoms (mainly O and N); molecules with TPSA < 140 A^2 typically show better membrane permeation.",
    },
    {
        "q": "Why must explicit hydrogen atoms be added before running 3D force-field geometry optimization (such as MMFF94)?",
        "options": [
            "Steric and electrostatic interactions depend on all atom positions including hydrogens",
            "RDKit cannot parse SMILES strings without explicit '[H]' symbols",
            "Hydrogens convert aromatic rings into linear alkynes",
            "Adding hydrogens reduces computation time to zero",
        ],
        "answer": "Steric and electrostatic interactions depend on all atom positions including hydrogens",
        "explanation": "Force fields like MMFF94 require explicit hydrogens to accurately compute bond lengths, angles, torsions, and van der Waals repulsions.",
    },
    {
        "q": "In SMILES notation, how are ring closures indicated (for example in Benzene, c1ccccc1)?",
        "options": [
            "Matching digits placed immediately after the atoms that close the ring",
            "Square brackets around the entire molecule",
            "Curly braces separating each carbon atom",
            "An asterisk (*) at both ends of the string",
        ],
        "answer": "Matching digits placed immediately after the atoms that close the ring",
        "explanation": "A digit (like '1') marks the opening atom of a ring and the matching digit marks the atom where the ring closes back.",
    },
]


# -----------------------------------------------------------------------------
# SESSION STATE & APPLE.COM UI/UX DESIGN SYSTEM
# -----------------------------------------------------------------------------

def init_session_state() -> None:
    """Initialize required session_state keys across all pages."""
    defaults: Dict[str, Any] = {
        "student_name": "Varun K",
        "reg_number": "RA2611026050122",
        "section": "AIML B",
        "visited_pages": {"Home"},
        "studied_molecules": {},
        "viz_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "prop_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "sim_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "batch_smiles": "\n".join(EXAMPLE_MOLECULES.values()),
        "quiz_index": 0,
        "quiz_answers": {},
        "quiz_submitted": {},
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def inject_apple_css() -> None:
    """Inject Apple.com minimalism CSS matching Mobbin reference specifications."""
    st.markdown(
        """
        <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap');

        :root {
            --apple-bg: #FBFBFD;
            --apple-surface: #FFFFFF;
            --apple-text-primary: #1D1D1F;
            --apple-text-secondary: #6E6E73;
            --apple-accent: #0071E3;
            --apple-accent-hover: #0077ED;
            --apple-border: #D2D2D7;
        }

        html, body, [class*="css"] {
            font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", "SF Pro Text", "Inter", sans-serif !important;
            color: #1D1D1F;
            -webkit-font-smoothing: antialiased;
        }

        .stApp {
            background-color: #FBFBFD;
            color: #1D1D1F;
        }

        .block-container {
            padding-top: 1.5rem !important;
            padding-bottom: 3rem !important;
            max-width: 1220px !important;
        }

        /* Sticky Translucent Top Navigation Bar (Apple Glassmorphism) */
        .apple-top-nav {
            position: sticky;
            top: 0;
            z-index: 999;
            background: rgba(251, 251, 253, 0.82);
            backdrop-filter: saturate(180%) blur(20px);
            -webkit-backdrop-filter: saturate(180%) blur(20px);
            border: 1px solid rgba(210, 210, 215, 0.75);
            border-radius: 16px;
            padding: 0.7rem 1.35rem;
            margin-bottom: 1.5rem;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 0.75rem;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.03);
        }

        .apple-brand {
            display: flex;
            align-items: center;
            gap: 0.55rem;
            font-weight: 700;
            font-size: 1.02rem;
            letter-spacing: -0.02em;
            color: #1D1D1F;
        }

        .apple-nav-pills {
            display: flex;
            align-items: center;
            gap: 0.4rem;
            flex-wrap: wrap;
        }

        .apple-nav-pill {
            font-size: 0.82rem;
            font-weight: 500;
            padding: 0.32rem 0.8rem;
            border-radius: 980px;
            color: #6E6E73;
            background: transparent;
            border: 1px solid transparent;
            transition: all 0.2s ease;
        }

        .apple-nav-pill.active {
            background: #1D1D1F;
            color: #FFFFFF;
            font-weight: 600;
        }

        /* Apple Hero Section */
        .apple-hero {
            background: linear-gradient(180deg, #FFFFFF 0%, #F5F5F7 100%);
            border: 1px solid #D2D2D7;
            border-radius: 20px;
            padding: 3.2rem 2.5rem;
            text-align: center;
            margin-bottom: 2rem;
            box-shadow: 0 12px 32px rgba(0, 0, 0, 0.04);
        }

        .apple-eyebrow {
            font-size: 0.82rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: #0071E3;
            margin-bottom: 0.65rem;
        }

        .apple-hero-title {
            font-size: 2.85rem;
            font-weight: 700;
            letter-spacing: -0.035em;
            line-height: 1.1;
            color: #1D1D1F;
            margin: 0 auto 0.9rem auto;
            max-width: 820px;
        }

        .apple-hero-subtitle {
            font-size: 1.15rem;
            font-weight: 400;
            line-height: 1.55;
            color: #6E6E73;
            max-width: 680px;
            margin: 0 auto 1.25rem auto;
        }

        /* Apple Minimalist Surface Cards */
        .apple-card {
            background-color: #FFFFFF;
            border: 1px solid #D2D2D7;
            border-radius: 18px;
            padding: 1.6rem;
            box-shadow: 0 4px 18px rgba(0, 0, 0, 0.03);
            transition: transform 0.25s cubic-bezier(0.16, 1, 0.3, 1), box-shadow 0.25s cubic-bezier(0.16, 1, 0.3, 1), border-color 0.25s ease;
            height: 100%;
            margin-bottom: 1rem;
        }

        .apple-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 14px 30px rgba(0, 0, 0, 0.06);
            border-color: #0071E3;
        }

        .apple-card-tag {
            display: inline-block;
            font-size: 0.75rem;
            font-weight: 600;
            color: #0071E3;
            background: rgba(0, 113, 227, 0.08);
            padding: 0.22rem 0.68rem;
            border-radius: 980px;
            margin-bottom: 0.65rem;
        }

        .apple-card h3 {
            font-size: 1.22rem;
            font-weight: 600;
            letter-spacing: -0.015em;
            color: #1D1D1F;
            margin: 0 0 0.45rem 0;
        }

        .apple-card p {
            font-size: 0.94rem;
            line-height: 1.55;
            color: #6E6E73;
            margin: 0;
        }

        /* Section Headlines */
        .apple-section-title {
            font-size: 1.85rem;
            font-weight: 700;
            letter-spacing: -0.025em;
            color: #1D1D1F;
            margin-top: 0.5rem;
            margin-bottom: 0.3rem;
        }

        .apple-section-lead {
            font-size: 1.02rem;
            color: #6E6E73;
            margin-bottom: 1.4rem;
        }

        /* Pill Buttons & Segmented Controls */
        .stButton > button, .stDownloadButton > button {
            border-radius: 980px !important;
            font-weight: 500 !important;
            padding: 0.5rem 1.25rem !important;
            border: 1px solid #D2D2D7 !important;
            transition: all 0.2s cubic-bezier(0.16, 1, 0.3, 1) !important;
        }

        .stButton > button[kind="primary"], .stDownloadButton > button[kind="primary"] {
            background-color: #0071E3 !important;
            color: #FFFFFF !important;
            border-color: #0071E3 !important;
        }

        .stButton > button:hover, .stDownloadButton > button:hover {
            transform: scale(1.015);
            box-shadow: 0 6px 16px rgba(0, 113, 227, 0.16);
        }

        /* Metric Cards */
        div[data-testid="stMetric"] {
            background-color: #FFFFFF;
            border: 1px solid #D2D2D7;
            border-radius: 16px;
            padding: 1rem 1.2rem;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.02);
        }

        /* Status Badges */
        .apple-badge-pass {
            display: inline-block;
            background-color: #E8F8EE;
            color: #116329;
            border: 1px solid #82E0A3;
            padding: 0.28rem 0.85rem;
            border-radius: 980px;
            font-weight: 600;
            font-size: 0.84rem;
        }

        .apple-badge-fail {
            display: inline-block;
            background-color: #FFEBEA;
            color: #C92A2A;
            border: 1px solid #FFA8A8;
            padding: 0.28rem 0.85rem;
            border-radius: 980px;
            font-weight: 600;
            font-size: 0.84rem;
        }

        /* Mandatory Apple Footer */
        .apple-footer {
            margin-top: 3.5rem;
            padding-top: 1.4rem;
            padding-bottom: 1.2rem;
            border-top: 1px solid #D2D2D7;
            text-align: center;
            color: #6E6E73;
            font-size: 0.85rem;
            font-weight: 400;
            letter-spacing: -0.005em;
        }

        @media (max-width: 768px) {
            .apple-hero-title { font-size: 2rem; }
            .apple-hero { padding: 2rem 1.2rem; }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_top_nav(active_page: str) -> None:
    """Render sticky translucent Apple-style top navigation bar and sidebar tracker."""
    init_session_state()
    st.session_state["visited_pages"].add(active_page)

    pages = [
        ("Home", "Home"),
        ("Molecule Viewer", "Viewer"),
        ("Property Calculator", "Properties"),
        ("Similarity Search", "Similarity"),
        ("Dataset Explorer", "Dataset"),
        ("About", "Theory & Quiz"),
    ]
    pills_html = "".join(
        f'<span class="apple-nav-pill {"active" if key == active_page else ""}">{label}</span>'
        for key, label in pages
    )

    st.markdown(
        f"""
        <div class="apple-top-nav">
            <div class="apple-brand">
                <span>🧪 Cheminformatics Virtual Lab</span>
                <span style="color:#6E6E73; font-weight:400; font-size:0.86rem;">| Apple Edition</span>
            </div>
            <div class="apple-nav-pills">
                {pills_html}
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.sidebar:
        st.markdown("### 🧪 Virtual Lab Profile")
        st.caption("SRM IST • Cheminformatics Suite")
        st.session_state["student_name"] = st.text_input(
            "Student Name",
            value=st.session_state["student_name"],
            placeholder="Varun K",
        )
        st.session_state["reg_number"] = st.text_input(
            "Registration Number",
            value=st.session_state["reg_number"],
            placeholder="RA2611026050122",
        )
        st.session_state["section"] = st.text_input(
            "Section",
            value=st.session_state["section"],
            placeholder="AIML B",
        )

        st.divider()
        visited_count = len(st.session_state["visited_pages"])
        total_modules = 6
        ratio = min(1.0, visited_count / total_modules)
        st.markdown("#### 📈 Lab Progress")
        st.progress(ratio, text=f"{visited_count} / {total_modules} Modules Explored ({ratio:.0%})")
        with st.expander("Completed Modules", expanded=False):
            for key, _ in pages:
                mark = "✅" if key in st.session_state["visited_pages"] else "⬜"
                st.write(f"{mark} {key}")

        with st.expander("🎯 Learning Objectives", expanded=False):
            for obj in LEARNING_OBJECTIVES:
                st.markdown(f"- {obj}")


def render_footer() -> None:
    """Render the mandatory centered Apple-style footer on every page."""
    st.markdown(
        f"""
        <div class="apple-footer">
            {FOOTER_TEXT}
        </div>
        """,
        unsafe_allow_html=True,
    )


# -----------------------------------------------------------------------------
# RDKit CHEMINFORMATICS FUNCTIONS (CACHED)
# -----------------------------------------------------------------------------

def validate_smiles(smiles: str) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate a SMILES string and return (Mol, friendly_error_message)."""
    cleaned = (smiles or "").strip()
    if not cleaned:
        return None, (
            "SMILES input is empty. Please enter a valid SMILES string such as "
            "`CCO` (Ethanol) or choose a preset molecule."
        )
    try:
        mol = Chem.MolFromSmiles(cleaned)
        if mol is None:
            return None, (
                f"Could not parse `{cleaned}` as a valid chemical structure. "
                "Check for unmatched parentheses `()`, unclosed ring digits (e.g. `c1ccccc1`), "
                "or invalid element symbols."
            )
        return mol, None
    except Exception as exc:
        return None, f"SMILES parsing error: {exc}. Try a canonical SMILES like `CC(=O)Oc1ccccc1C(=O)O`."


def _draw_2d_fallback_pillow(mol: Chem.Mol, width: int = 500, height: int = 380) -> bytes:
    """Pure Pillow 2D structure renderer fallback when X11 libraries are unavailable."""
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)

    if n_atoms == 0:
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()

    coords = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(n_atoms)]
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    span_x = max(max_x - min_x, 1.5)
    span_y = max(max_y - min_y, 1.5)
    scale = min((width - 100) / span_x, (height - 90) / span_y, 68.0)
    cx = width / 2.0 - ((min_x + max_x) / 2.0) * scale
    cy = height / 2.0 + ((min_y + max_y) / 2.0) * scale

    def to_px(idx: int) -> Tuple[float, float]:
        x, y = coords[idx]
        return (cx + x * scale, cy - y * scale)

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        p1 = to_px(i)
        p2 = to_px(j)
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        dist = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / dist, dx / dist
        btype = bond.GetBondTypeAsDouble()

        if btype >= 2.0 or bond.GetIsAromatic():
            offset = 3.5
            draw.line([(p1[0] + nx * offset, p1[1] + ny * offset), (p2[0] + nx * offset, p2[1] + ny * offset)], fill="#1D1D1F", width=3)
            draw.line([(p1[0] - nx * offset, p1[1] - ny * offset), (p2[0] - nx * offset, p2[1] - ny * offset)], fill="#1D1D1F", width=3)
        else:
            draw.line([p1, p2], fill="#1D1D1F", width=3)

    atom_colors = {"O": "#E03131", "N": "#0071E3", "S": "#D97706", "F": "#0CA678", "Cl": "#2F9E44", "Br": "#7950F2"}
    font = ImageFont.load_default()
    for idx in range(n_atoms):
        atom = mol.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        if sym != "C" or n_atoms <= 2:
            px_x, px_y = to_px(idx)
            r = 12
            draw.ellipse([px_x - r, px_y - r, px_x + r, px_y + r], fill="white")
            col = atom_colors.get(sym, "#1D1D1F")
            draw.text((px_x - 4 * len(sym), px_y - 6), sym, fill=col, font=font)

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


@st.cache_data(show_spinner=False)
def render_2d_image(smiles: str, width: int = 500, height: int = 380) -> Optional[bytes]:
    """Generate a high-clarity PNG byte buffer of the 2D chemical structure."""
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    if HAS_RDKIT_DRAW and RDKitDraw is not None:
        try:
            img: Image.Image = RDKitDraw.MolToImage(mol, size=(width, height))
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            return buf.getvalue()
        except Exception:
            pass
    return _draw_2d_fallback_pillow(mol, width=width, height=height)


@st.cache_data(show_spinner=False)
def build_3d_molblock(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Protonate and embed a 3D conformer optimized via MMFF94/UFF."""
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None, "Invalid SMILES for 3D structure generation."
    try:
        mol_3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        res = AllChem.EmbedMolecule(mol_3d, params)
        if res == -1:
            res = AllChem.EmbedMolecule(mol_3d, useRandomCoords=True, randomSeed=42)
        if res != -1:
            try:
                AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=200)
            except Exception:
                AllChem.UFFOptimizeMolecule(mol_3d, maxIters=200)
        return Chem.MolToMolBlock(mol_3d), None
    except Exception as exc:
        return None, f"3D conformer embedding failed ({exc}). Try a smaller organic molecule."


def render_3d_viewer(view: Any, height: int = 440, width: int = 720) -> None:
    """Render py3Dmol viewer via stmol or Streamlit HTML component fallback."""
    if HAS_STMOL and _stmol_showmol is not None:
        try:
            _stmol_showmol(view, height=height, width=width)
            return
        except Exception:
            pass
    components.html(view._make_html(), height=height + 20, width=width)


@st.cache_data(show_spinner=False)
def calculate_molecular_properties(smiles: str, custom_name: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Compute physicochemical descriptors and Lipinski Rule of Five status."""
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None

    canonical = Chem.MolToSmiles(mol)
    mw = round(Descriptors.MolWt(mol), 2)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    logp = round(Crippen.MolLogP(mol), 2)
    tpsa = round(Descriptors.TPSA(mol), 2)
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rot_bonds = int(Lipinski.NumRotatableBonds(mol))
    rings = int(Lipinski.RingCount(mol))
    heavy_atoms = int(mol.GetNumHeavyAtoms())
    fsp3 = round(rdMolDescriptors.CalcFractionCSP3(mol), 2)

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    lipinski_pass = violations <= 1

    name_lookup = {
        Chem.MolToSmiles(Chem.MolFromSmiles(item["SMILES"])): item["Name"]
        for item in BENCHMARK_LIBRARY
        if Chem.MolFromSmiles(item["SMILES"]) is not None
    }
    mol_name = custom_name or name_lookup.get(canonical, canonical)

    return {
        "Name": mol_name,
        "SMILES": canonical,
        "Formula": formula,
        "Molecular Weight": mw,
        "LogP": logp,
        "TPSA": tpsa,
        "H-bond Donors": hbd,
        "H-bond Acceptors": hba,
        "Rotatable Bonds": rot_bonds,
        "Ring Count": rings,
        "Heavy Atoms": heavy_atoms,
        "Fraction Csp3": fsp3,
        "Lipinski Violations": violations,
        "Lipinski Check": "Pass" if lipinski_pass else "Fail",
    }


def record_studied_molecule(props: Dict[str, Any]) -> None:
    """Save calculated molecule properties in st.session_state for session portfolio."""
    init_session_state()
    st.session_state["studied_molecules"][props["SMILES"]] = props


@st.cache_data(show_spinner=False)
def compute_similarity_search(
    query_smiles: str,
    fp_type: str = "Morgan (ECFP4)",
    radius: int = 2,
    n_bits: int = 2048,
    custom_library: Optional[List[Dict[str, str]]] = None,
) -> pd.DataFrame:
    """Compute Tanimoto similarity between a query SMILES and a library of molecules."""
    query_mol = Chem.MolFromSmiles(query_smiles.strip())
    if query_mol is None:
        return pd.DataFrame()

    def _get_fp(m: Chem.Mol):
        if fp_type.startswith("MACCS"):
            return rdMolDescriptors.GetMACCSKeysFingerprint(m)
        if fp_type.startswith("RDKit"):
            return Chem.RDKFingerprint(m, fpSize=n_bits)
        return AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=n_bits)

    q_fp = _get_fp(query_mol)
    lib = custom_library if custom_library is not None else BENCHMARK_LIBRARY

    rows: List[Dict[str, Any]] = []
    for entry in lib:
        smi = entry.get("SMILES", "").strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        t_fp = _get_fp(mol)
        sim = float(DataStructs.TanimotoSimilarity(q_fp, t_fp))
        on_bits_q = int(q_fp.GetNumOnBits())
        on_bits_t = int(t_fp.GetNumOnBits())
        props = calculate_molecular_properties(smi, custom_name=entry.get("Name"))
        if props:
            rows.append(
                {
                    "Name": props["Name"],
                    "Category": entry.get("Category", "Library Compound"),
                    "SMILES": props["SMILES"],
                    "Formula": props["Formula"],
                    "Tanimoto Similarity": round(sim, 4),
                    "Molecular Weight": props["Molecular Weight"],
                    "LogP": props["LogP"],
                    "TPSA": props["TPSA"],
                    "Active Bits": on_bits_t,
                    "Lipinski Check": props["Lipinski Check"],
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(by="Tanimoto Similarity", ascending=False).reset_index(drop=True)
    return df


def generate_pdf_report(
    title: str,
    student_name: str,
    reg_number: str,
    properties_list: List[Dict[str, Any]],
    assessment_score: Optional[Tuple[int, int]] = None,
) -> bytes:
    """Build an Apple-inspired academic PDF report with SRM IST Division of Chemistry header and Dr. P Pachamuthu faculty attribution."""
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=letter,
        rightMargin=40,
        leftMargin=40,
        topMargin=36,
        bottomMargin=36,
    )
    styles = getSampleStyleSheet()

    inst_style = ParagraphStyle(
        "SRMInstHeader",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=12,
        leading=15,
        alignment=1,  # Center
        textColor=colors.HexColor("#1D1D1F"),
        spaceAfter=2,
    )
    dept_style = ParagraphStyle(
        "SRMDeptHeader",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9.5,
        leading=13,
        alignment=1,  # Center
        textColor=colors.HexColor("#6E6E73"),
        spaceAfter=10,
    )
    title_style = ParagraphStyle(
        "AppleTitle",
        parent=styles["Heading1"],
        fontName="Helvetica-Bold",
        fontSize=16,
        leading=20,
        alignment=1,  # Center
        textColor=colors.HexColor("#0071E3"),
        spaceBefore=4,
        spaceAfter=12,
    )
    subtitle_style = ParagraphStyle(
        "AppleSub",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=12,
        leading=16,
        textColor=colors.HexColor("#1D1D1F"),
        spaceBefore=12,
        spaceAfter=6,
    )
    meta_style = ParagraphStyle(
        "AppleMeta",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9.5,
        leading=14,
        textColor=colors.HexColor("#1D1D1F"),
        spaceAfter=6,
    )
    faculty_style = ParagraphStyle(
        "FacultyBox",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9,
        leading=13,
        textColor=colors.HexColor("#1D1D1F"),
        spaceBefore=6,
        spaceAfter=8,
    )
    footer_pdf_style = ParagraphStyle(
        "FooterPDF",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=8.5,
        leading=12,
        alignment=1,
        textColor=colors.HexColor("#6E6E73"),
        spaceBefore=14,
    )

    story: List[Any] = []

    # 1. SRM IST Institutional Header
    story.append(Paragraph("SRM INSTITUTE OF SCIENCE AND TECHNOLOGY, TIRUCHIRAPPALLI", inst_style))
    story.append(
        Paragraph(
            "Faculty of Engineering and Technology &nbsp;|&nbsp; School of Sciences &nbsp;|&nbsp; Division of Chemistry",
            dept_style,
        )
    )
    story.append(Paragraph(f"Cheminformatics Virtual Laboratory — {title}", title_style))

    # 2. Student & Session Metadata Table
    section_val = st.session_state.get("section", "AIML B") if hasattr(st, "session_state") else "AIML B"
    story.append(
        Paragraph(
            f"<b>Student Name:</b> {student_name or 'Varun K'} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Registration No:</b> {reg_number or 'RA2611026050122'} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Section:</b> {section_val} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Date:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            meta_style,
        )
    )
    story.append(Spacer(1, 10))

    # 3. Molecular Descriptors & Lipinski Summary Table
    if properties_list:
        story.append(Paragraph("Molecular Descriptors & Lipinski Rule-of-Five Summary", subtitle_style))
        headers = ["Name / SMILES", "Formula", "MW", "LogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "Lipinski"]
        table_data = [headers]
        for item in properties_list:
            label = str(item.get("Name", item.get("SMILES", "")))
            if len(label) > 18:
                label = label[:15] + "..."
            table_data.append(
                [
                    label,
                    str(item.get("Formula", "")),
                    f"{float(item.get('Molecular Weight', 0)):.2f}",
                    f"{float(item.get('LogP', 0)):.2f}",
                    f"{float(item.get('TPSA', 0)):.1f}",
                    str(item.get("H-bond Donors", 0)),
                    str(item.get("H-bond Acceptors", 0)),
                    str(item.get("Rotatable Bonds", 0)),
                    str(item.get("Ring Count", 0)),
                    str(item.get("Lipinski Check", "")),
                ]
            )
        t = Table(table_data, repeatRows=1)
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1D1D1F")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                    ("ALIGN", (1, 0), (-1, -1), "CENTER"),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#D2D2D7")),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#FBFBFD")]),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 6),
                ]
            )
        )
        story.append(t)
        story.append(Spacer(1, 12))

    # 4. Assessment Score Section (if applicable)
    if assessment_score is not None:
        score, total = assessment_score
        pct = (score / total * 100) if total else 0.0
        status = "PASS (Distinction)" if pct >= 80 else ("PASS" if pct >= 60 else "NEEDS REVIEW")
        story.append(Paragraph("Self-Assessment Quiz Performance", subtitle_style))
        story.append(
            Paragraph(
                f"<b>Quiz Score:</b> {score} / {total} ({pct:.0f}%) &nbsp;&nbsp;|&nbsp;&nbsp; "
                f"<b>Achievement Status:</b> {status}",
                meta_style,
            )
        )
        story.append(Spacer(1, 10))

    # 5. Faculty Attribution Block (Dr. P Pachamuthu)
    story.append(Paragraph("Faculty Attribution & Course Coordinator", subtitle_style))
    story.append(
        Paragraph(
            "<b>Dr. P Pachamuthu</b> — Professor, Division of Chemistry, School of Sciences, "
            "Faculty of Engineering and Technology, SRM Institute of Science and Technology, Tiruchirappalli "
            "(Contact: pachamuthup@srmist.edu.in)",
            faculty_style,
        )
    )

    # 6. Mandatory Student Footer
    story.append(Paragraph(FOOTER_TEXT, footer_pdf_style))
    doc.build(story)
    return buffer.getvalue()

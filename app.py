"""
Cheminformatics Virtual Laboratory
Interactive Streamlit multi-page web app for undergraduate chemistry students.
Python 3.11 | RDKit | py3Dmol | stmol | Plotly | ReportLab
"""

from __future__ import annotations

import io
import math
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import plotly.express as px
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

# Safe import for RDKit 2D Draw (requires libXrender1 on Linux containers)
try:
    from rdkit.Chem import Draw as RDKitDraw
    HAS_RDKIT_DRAW = True
except Exception:
    RDKitDraw = None
    HAS_RDKIT_DRAW = False

# Safe import for stmol with native Streamlit HTML component fallback
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
from streamlit_option_menu import option_menu

# -----------------------------------------------------------------------------
# CONSTANTS & REFERENCE DATA
# -----------------------------------------------------------------------------

PAGES: List[str] = [
    "Home",
    "Theory",
    "Molecular Visualization",
    "Molecular Properties Report",
    "Structure-Property Analysis",
    "Assessment",
    "Final Project Report",
]

PAGE_ICONS: List[str] = [
    "house-door-fill",
    "book-half",
    "eyedropper",
    "clipboard2-data-fill",
    "bar-chart-line-fill",
    "patch-question-fill",
    "file-earmark-pdf-fill",
]

LEARNING_OBJECTIVES: List[str] = [
    "SMILES notation",
    "Molecular representation",
    "2D visualization",
    "3D visualization",
    "Molecular descriptors",
    "Structure-property relationships",
]

EXAMPLE_MOLECULES: Dict[str, str] = {
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Phenol": "Oc1ccccc1",
    "Acetic Acid": "CC(=O)O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
}

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
        "options": [
            "<= 500 Da",
            "<= 100 Da",
            ">= 500 Da",
            "Exactly 250 Da",
        ],
        "answer": "<= 500 Da",
        "explanation": "Lipinski's criteria stipulate MW <= 500 Da, LogP <= 5, H-bond donors <= 5, and H-bond acceptors <= 10 for good oral bioavailability.",
    },
    {
        "q": "What is the maximum number of Hydrogen Bond Donors (HBD) allowed under Lipinski's Rule of Five?",
        "options": ["5", "10", "2", "15"],
        "answer": "5",
        "explanation": "Lipinski's rule specifies no more than 5 hydrogen bond donors (total N-H and O-H bonds).",
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
# SESSION STATE & THEME STYLING
# -----------------------------------------------------------------------------

def init_session_state() -> None:
    """Initialize required session_state keys across all lab pages."""
    defaults: Dict[str, Any] = {
        "current_page": "Home",
        "student_name": "",
        "reg_number": "",
        "dark_mode": False,
        "completed_pages": {"Home"},
        "studied_molecules": {},
        "viz_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "prop_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "batch_smiles": "\n".join(EXAMPLE_MOLECULES.values()),
        "quiz_index": 0,
        "quiz_answers": {},
        "quiz_submitted": {},
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def inject_custom_css(dark_mode: bool = False) -> None:
    """Apply Inter font, custom academic palette, rounded cards, and badges."""
    bg_color = "#0F172A" if dark_mode else "#F8FAFC"
    card_bg = "#1E293B" if dark_mode else "#FFFFFF"
    text_color = "#F8FAFC" if dark_mode else "#0F172A"
    muted_text = "#94A3B8" if dark_mode else "#475569"
    border_color = "#334155" if dark_mode else "#E2E8F0"

    st.markdown(
        f"""
        <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

        html, body, [class*="css"] {{
            font-family: 'Inter', sans-serif;
            color: {text_color};
        }}
        .stApp {{
            background-color: {bg_color};
            color: {text_color};
        }}
        .hero-card {{
            background: linear-gradient(135deg, #1E3A8A 0%, #0F296B 65%, #14B8A6 100%);
            color: #FFFFFF;
            padding: 2.2rem 2rem;
            border-radius: 16px;
            box-shadow: 0 10px 25px -5px rgba(30, 58, 138, 0.25);
            margin-bottom: 1.5rem;
        }}
        .hero-title {{
            font-size: 2.1rem;
            font-weight: 700;
            margin: 0.3rem 0;
            letter-spacing: -0.02em;
        }}
        .hero-tagline {{
            font-size: 1.05rem;
            opacity: 0.92;
            margin-bottom: 0.5rem;
        }}
        .lab-card {{
            background-color: {card_bg};
            border: 1px solid {border_color};
            border-radius: 14px;
            padding: 1.35rem;
            box-shadow: 0 4px 12px rgba(15, 23, 42, 0.05);
            margin-bottom: 1rem;
            height: 100%;
        }}
        .lab-card h4 {{
            color: {"#38BDF8" if dark_mode else "#1E3A8A"};
            margin-top: 0;
            margin-bottom: 0.5rem;
            font-weight: 600;
        }}
        .lab-card p {{
            color: {muted_text};
            font-size: 0.93rem;
            line-height: 1.5;
            margin-bottom: 0;
        }}
        .badge-pass {{
            display: inline-block;
            background-color: #CCFBF1;
            color: #0F766E;
            border: 1px solid #14B8A6;
            padding: 0.3rem 0.85rem;
            border-radius: 9999px;
            font-weight: 600;
            font-size: 0.88rem;
        }}
        .badge-fail {{
            display: inline-block;
            background-color: #FEE2E2;
            color: #B91C1C;
            border: 1px solid #EF4444;
            padding: 0.3rem 0.85rem;
            border-radius: 9999px;
            font-weight: 600;
            font-size: 0.88rem;
        }}
        .footer-box {{
            text-align: center;
            padding: 1.25rem;
            margin-top: 2rem;
            border-top: 1px solid {border_color};
            color: {muted_text};
            font-size: 0.88rem;
        }}
        @media (max-width: 768px) {{
            .hero-title {{ font-size: 1.55rem; }}
            .hero-card {{ padding: 1.4rem 1.2rem; }}
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )


# -----------------------------------------------------------------------------
# CHEMISTRY HELPER FUNCTIONS (Cached & Error-Safe)
# -----------------------------------------------------------------------------

def validate_smiles(smiles: str) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """Validate a SMILES string and return (Mol, friendly_error_message)."""
    cleaned = (smiles or "").strip()
    if not cleaned:
        return None, (
            "SMILES input is empty. Suggested fix: Enter a valid SMILES string "
            "such as `CCO` (Ethanol) or click one of the example buttons."
        )
    try:
        mol = Chem.MolFromSmiles(cleaned)
        if mol is None:
            return None, (
                f"Could not parse `{cleaned}` as a valid chemical structure. "
                "Suggested fix: Check for unmatched parentheses `()`, unclosed ring numbers "
                "(e.g., use `c1ccccc1` with matching `1`s), or invalid atomic symbols."
            )
        return mol, None
    except Exception as exc:
        return None, (
            f"Unexpected SMILES parsing error: {exc}. "
            "Suggested fix: Try a simpler canonical SMILES like `CC(=O)O`."
        )


def _draw_2d_fallback_pillow(mol: Chem.Mol, width: int = 450, height: int = 340) -> bytes:
    """Pure Pillow 2D structure renderer using RDKit 2D coordinates (no libXrender needed)."""
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
    scale = min((width - 100) / span_x, (height - 90) / span_y, 65.0)
    cx = width / 2.0 - ((min_x + max_x) / 2.0) * scale
    cy = height / 2.0 + ((min_y + max_y) / 2.0) * scale

    def to_px(idx: int) -> Tuple[float, float]:
        x, y = coords[idx]
        return (cx + x * scale, cy - y * scale)

    # Draw bonds
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
            draw.line([(p1[0] + nx * offset, p1[1] + ny * offset), (p2[0] + nx * offset, p2[1] + ny * offset)], fill="#1E293B", width=3)
            draw.line([(p1[0] - nx * offset, p1[1] - ny * offset), (p2[0] - nx * offset, p2[1] - ny * offset)], fill="#1E293B", width=3)
        else:
            draw.line([p1, p2], fill="#1E293B", width=3)

    # Draw heteroatom labels
    atom_colors = {"O": "#DC2626", "N": "#2563EB", "S": "#D97706", "F": "#0D9488", "Cl": "#16A34A", "Br": "#9333EA"}
    font = ImageFont.load_default()
    for idx in range(n_atoms):
        atom = mol.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        if sym != "C" or n_atoms <= 2:
            px_x, px_y = to_px(idx)
            r = 12
            draw.ellipse([px_x - r, px_y - r, px_x + r, px_y + r], fill="white")
            col = atom_colors.get(sym, "#0F172A")
            draw.text((px_x - 4 * len(sym), px_y - 6), sym, fill=col, font=font)

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


@st.cache_data(show_spinner=False)
def render_2d_image(smiles: str, width: int = 450, height: int = 340) -> Optional[bytes]:
    """Generate a PNG byte buffer of the 2D molecular depiction with X11-safe fallback."""
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
    """Add explicit hydrogens and optimize 3D geometry via MMFF94/UFF."""
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None, "Invalid SMILES for 3D generation."
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
        return None, (
            f"3D conformer generation failed ({exc}). "
            "Suggested fix: Try a smaller organic molecule without disconnected salts."
        )


def render_3d_viewer(view: Any, height: int = 450, width: int = 700) -> None:
    """Render py3Dmol view using stmol or native Streamlit HTML component fallback."""
    if HAS_STMOL and _stmol_showmol is not None:
        try:
            _stmol_showmol(view, height=height, width=width)
            return
        except Exception:
            pass
    components.html(view._make_html(), height=height + 20, width=width)


@st.cache_data(show_spinner=False)
def calculate_molecular_properties(smiles: str) -> Optional[Dict[str, Any]]:
    """Compute core RDKit descriptors and Lipinski Rule of Five compliance."""
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

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    lipinski_pass = violations <= 1

    name_lookup = {Chem.MolToSmiles(Chem.MolFromSmiles(v)): k for k, v in EXAMPLE_MOLECULES.items()}
    mol_name = name_lookup.get(canonical, canonical)

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
        "Lipinski Violations": violations,
        "Lipinski Check": "Pass" if lipinski_pass else "Fail",
    }


def record_studied_molecule(props: Dict[str, Any]) -> None:
    """Save calculated molecule properties in st.session_state for the final report."""
    st.session_state["studied_molecules"][props["SMILES"]] = props


def generate_pdf_report(
    title: str,
    student_name: str,
    reg_number: str,
    properties_list: List[Dict[str, Any]],
    assessment_score: Optional[Tuple[int, int]] = None,
) -> bytes:
    """Build a clean academic PDF report using ReportLab."""
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=letter,
        rightMargin=40,
        leftMargin=40,
        topMargin=40,
        bottomMargin=40,
    )
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "LabTitle",
        parent=styles["Heading1"],
        fontSize=18,
        textColor=colors.HexColor("#1E3A8A"),
        spaceAfter=10,
    )
    subtitle_style = ParagraphStyle(
        "LabSub",
        parent=styles["Heading2"],
        fontSize=13,
        textColor=colors.HexColor("#14B8A6"),
        spaceBefore=12,
        spaceAfter=6,
    )
    normal_style = styles["Normal"]

    story: List[Any] = []
    story.append(Paragraph(f"Cheminformatics Virtual Laboratory — {title}", title_style))
    story.append(
        Paragraph(
            f"<b>Student Name:</b> {student_name or 'Not Provided'} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Registration No:</b> {reg_number or 'Not Provided'} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Date:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            normal_style,
        )
    )
    story.append(Spacer(1, 14))

    if properties_list:
        story.append(Paragraph("Molecular Descriptors & Lipinski Summary", subtitle_style))
        headers = ["Name / SMILES", "Formula", "MW", "LogP", "TPSA", "HBD", "HBA", "RotB", "Rings", "Lipinski"]
        table_data = [headers]
        for item in properties_list:
            label = item.get("Name", item["SMILES"])
            if len(label) > 18:
                label = label[:15] + "..."
            table_data.append(
                [
                    label,
                    str(item["Formula"]),
                    f"{item['Molecular Weight']:.2f}",
                    f"{item['LogP']:.2f}",
                    f"{item['TPSA']:.1f}",
                    str(item["H-bond Donors"]),
                    str(item["H-bond Acceptors"]),
                    str(item["Rotatable Bonds"]),
                    str(item["Ring Count"]),
                    str(item["Lipinski Check"]),
                ]
            )
        t = Table(table_data, repeatRows=1)
        t.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1E3A8A")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                    ("ALIGN", (1, 0), (-1, -1), "CENTER"),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#CBD5E1")),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8FAFC")]),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 6),
                ]
            )
        )
        story.append(t)
        story.append(Spacer(1, 14))

    if assessment_score is not None:
        score, total = assessment_score
        pct = (score / total * 100) if total else 0.0
        status = "PASS (Distinction)" if pct >= 80 else ("PASS" if pct >= 50 else "NEEDS REVIEW")
        story.append(Paragraph("Assessment Performance", subtitle_style))
        story.append(
            Paragraph(
                f"<b>Quiz Score:</b> {score} / {total} ({pct:.0f}%) &nbsp;&nbsp;|&nbsp;&nbsp; "
                f"<b>Achievement Status:</b> {status}",
                normal_style,
            )
        )

    doc.build(story)
    return buffer.getvalue()


# -----------------------------------------------------------------------------
# SIDEBAR NAVIGATION & STUDENT TRACKER
# -----------------------------------------------------------------------------

def render_sidebar() -> str:
    """Render sidebar with option_menu, student info, progress bar, and objectives."""
    with st.sidebar:
        st.markdown("### 🧪 Cheminformatics Lab")
        st.caption("Undergraduate Virtual Laboratory")

        default_idx = PAGES.index(st.session_state["current_page"])
        selected_page = option_menu(
            menu_title="Lab Navigation",
            options=PAGES,
            icons=PAGE_ICONS,
            menu_icon="mortarboard-fill",
            default_index=default_idx,
            styles={
                "container": {"padding": "0.3rem", "background-color": "transparent"},
                "icon": {"color": "#14B8A6", "font-size": "16px"},
                "nav-link": {
                    "font-size": "14px",
                    "text-align": "left",
                    "margin": "2px 0",
                    "border-radius": "8px",
                },
                "nav-link-selected": {
                    "background-color": "#1E3A8A",
                    "color": "#FFFFFF",
                    "font-weight": "600",
                },
            },
        )

        if selected_page != st.session_state["current_page"]:
            st.session_state["current_page"] = selected_page
            st.session_state["completed_pages"].add(selected_page)
            st.rerun()

        st.session_state["completed_pages"].add(selected_page)

        st.divider()
        st.markdown("#### 👤 Student Information")
        st.session_state["student_name"] = st.text_input(
            "Student Name",
            value=st.session_state["student_name"],
            placeholder="e.g., Priya Sharma",
            help="Enter your full name to personalize all downloadable PDF and CSV lab reports.",
        )
        st.session_state["reg_number"] = st.text_input(
            "Registration Number",
            value=st.session_state["reg_number"],
            placeholder="e.g., CHM2026-104",
            help="Enter your university roll or registration number for lab grading.",
        )

        st.divider()
        st.markdown("#### 📈 Lab Progress Tracker")
        done_count = len(st.session_state["completed_pages"])
        total_pages = len(PAGES)
        progress_ratio = done_count / total_pages
        st.progress(progress_ratio, text=f"{done_count} / {total_pages} Modules Visited ({progress_ratio:.0%})")
        with st.expander("View Completed Modules"):
            for page in PAGES:
                mark = "✅" if page in st.session_state["completed_pages"] else "⬜"
                st.write(f"{mark} {page}")

        st.divider()
        with st.expander("🎯 Learning Objectives", expanded=False):
            for obj in LEARNING_OBJECTIVES:
                st.markdown(f"- {obj}")

        st.divider()
        st.session_state["dark_mode"] = st.toggle(
            "🌙 Dark Mode",
            value=st.session_state["dark_mode"],
            help="Toggle between clean academic light mode and high-contrast dark mode.",
        )

    return st.session_state["current_page"]


# -----------------------------------------------------------------------------
# PAGE 1: HOME
# -----------------------------------------------------------------------------

def render_home() -> None:
    """Render the Home page with hero banner, feature cards, example table, and footer."""
    st.markdown(
        """
        <div class="hero-card">
            <div style="font-size: 0.9rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.08em; color: #5EEAD4;">
                🏛️ Department of Chemistry • Virtual Laboratory
            </div>
            <div class="hero-title">Cheminformatics Virtual Laboratory</div>
            <div class="hero-tagline">
                An interactive, modern, beginner-friendly workspace for exploring molecular representations,
                2D/3D conformations, physicochemical descriptors, and drug-likeness rules.
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    col_cta, col_info = st.columns([1, 3])
    with col_cta:
        if st.button(
            "🚀 Start Learning",
            type="primary",
            use_container_width=True,
            help="Jump directly to the Theory module to begin the guided cheminformatics workflow.",
        ):
            st.session_state["current_page"] = "Theory"
            st.session_state["completed_pages"].add("Theory")
            st.toast("Welcome to the Theory module!", icon="📘")
            st.rerun()
    with col_info:
        st.info(
            "💡 **Tip:** Fill in your **Student Name** and **Registration Number** in the left sidebar "
            "so your downloadable lab reports automatically include your details."
        )

    st.subheader("Welcome to Computational Chemistry & Cheminformatics")
    st.write(
        "Cheminformatics bridges experimental chemistry and data science. In this virtual laboratory, "
        "you will learn how computers encode chemical structures as text strings, visualize molecules "
        "in two and three dimensions, compute quantitative descriptors that govern biological activity, "
        "and evaluate drug-likeness using industry-standard rules."
    )

    st.subheader("What Can Students Do?")
    features = [
        (
            "📘 Master Core Theory",
            "Explore concise, interactive explanations of SMILES notation, molecular graphs, topological descriptors, and Lipinski's Rule of Five.",
        ),
        (
            "🔬 2D & 3D Visualization",
            "Convert SMILES strings into crisp RDKit 2D chemical drawings and interactive, force-field optimized 3D models via py3Dmol.",
        ),
        (
            "📊 Molecular Properties Report",
            "Compute MW, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, and export formatted CSV and PDF lab sheets.",
        ),
        (
            "📈 Structure-Property Analysis",
            "Compare multiple molecules side-by-side with interactive Plotly bar charts and MW vs. LogP chemical space scatter plots.",
        ),
        (
            "✅ Self-Assessment Quiz",
            "Test your understanding with a 10-question interactive quiz featuring instant explanations and achievement badges.",
        ),
        (
            "🎓 Final Project Portfolio",
            "Compile all molecules studied during your session along with your quiz score into a comprehensive PDF project report.",
        ),
    ]

    for row_start in range(0, len(features), 3):
        cols = st.columns(3)
        for col, (f_title, f_desc) in zip(cols, features[row_start : row_start + 3]):
            with col:
                st.markdown(
                    f"""
                    <div class="lab-card">
                        <h4>{f_title}</h4>
                        <p>{f_desc}</p>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

    st.subheader("Reference Example Molecules")
    st.caption("Use these benchmark SMILES strings across the Visualization and Properties modules.")
    example_rows = []
    for name, smi in EXAMPLE_MOLECULES.items():
        props = calculate_molecular_properties(smi)
        if props:
            example_rows.append(
                {
                    "Molecule Name": name,
                    "SMILES Notation": smi,
                    "Molecular Formula": props["Formula"],
                    "Molecular Weight (g/mol)": props["Molecular Weight"],
                    "LogP": props["LogP"],
                    "Lipinski Rule": props["Lipinski Check"],
                }
            )
    st.dataframe(pd.DataFrame(example_rows), use_container_width=True, hide_index=True)

    st.markdown(
        """
        <div class="footer-box">
            <strong>Developed for Undergraduate Chemistry Education</strong> •
            Powered by Python 3.11, Streamlit, RDKit, py3Dmol &amp; Plotly
        </div>
        """,
        unsafe_allow_html=True,
    )


# -----------------------------------------------------------------------------
# PAGE 2: THEORY
# -----------------------------------------------------------------------------

def render_theory() -> None:
    """Render the Theory page using expanders with short text and examples."""
    st.header("📘 Cheminformatics Theory & Foundations")
    st.write(
        "Click each topic below to review the theoretical principles before conducting "
        "your virtual experiments."
    )

    with st.expander("1. SMILES (Simplified Molecular Input Line Entry System)", expanded=True):
        st.markdown(
            """
            **Concept:** SMILES is a line notation that represents chemical structures using short ASCII strings
            so computers can rapidly store, search, and analyze millions of molecules.

            **Key Rules:**
            - **Atoms:** Represented by atomic symbols (`C`, `N`, `O`, `S`, `F`, `Cl`, `Br`). Hydrogens are implicit.
            - **Bonds:** Single bonds are implicit (`CC` = ethane), double bonds use `=`, triple bonds use `#`.
            - **Branches:** Enclosed in parentheses, e.g., `CC(=O)O` for acetic acid where `=O` branches off the second carbon.
            - **Rings & Aromaticity:** Matching digits close rings, and lowercase letters indicate aromatic atoms (`c1ccccc1` = benzene).

            **Worked Examples:**
            - **Ethanol (`CCO`):** Two aliphatic carbons connected to a hydroxyl oxygen.
            - **Phenol (`Oc1ccccc1`):** A hydroxyl group attached to a 6-membered aromatic benzene ring.
            """
        )

    with st.expander("2. Molecular Graphs", expanded=False):
        st.markdown(
            """
            **Concept:** Inside cheminformatics toolkits like RDKit, every molecule is stored mathematically as a
            **labeled graph** $G = (V, E)$:
            - **Vertices / Nodes ($V$):** Represent **atoms** storing properties like atomic number, formal charge, hybridization, and chirality.
            - **Edges ($E$):** Represent **chemical bonds** storing bond order (single, double, triple, aromatic) and stereochemistry ($E/Z$).

            **Example:**
            - **Acetic Acid (`CC(=O)O`)** has 4 heavy-atom nodes (2 Carbons, 2 Oxygens) connected by 3 edges (one C–C single bond, one C=O double bond, and one C–O single bond).
            """
        )

    with st.expander("3. Molecular Descriptors", expanded=False):
        st.markdown(
            """
            **Concept:** Molecular descriptors are numerical values extracted from a chemical structure to quantify
            its physical, chemical, or topological traits:
            - **Molecular Weight (MW):** Sum of atomic weights of all atoms in the molecule (g/mol or Da).
            - **LogP (Octanol-Water Partition Coefficient):** Measures lipophilicity ($\log_{10}([\\text{Octanol}]/[\\text{Water}])$). Higher LogP means more hydrophobic.
            - **TPSA (Topological Polar Surface Area):** Surface area ($\text{\AA}^2$) occupied by polar nitrogen and oxygen atoms and their attached hydrogens.
            - **Rotatable Bonds:** Non-ring single bonds bound to non-terminal heavy atoms; measures conformational flexibility.

            **Example:**
            - **Aspirin (`CC(=O)Oc1ccccc1C(=O)O`)** has $\\text{MW} = 180.16\\text{ g/mol}$, $\\text{LogP} = 1.31$, and $\\text{TPSA} = 63.60\\text{ \AA}^2$.
            """
        )

    with st.expander("4. Lipinski's Rule of Five", expanded=False):
        st.markdown(
            """
            **Concept:** Proposed by Christopher Lipinski (1997), the **Rule of Five (Ro5)** evaluates whether a compound
            has physicochemical properties favorable for **oral bioavailability** in humans.

            **The Four Thresholds (multiples of 5):**
            1. **Molecular Weight (MW)** $\\le 500\\text{ Da}$
            2. **Lipophilicity (LogP)** $\\le 5$
            3. **Hydrogen Bond Donors (HBD)** $\\le 5$ (sum of $\\text{O-H}$ and $\\text{N-H}$ groups)
            4. **Hydrogen Bond Acceptors (HBA)** $\\le 10$ (sum of $\\text{N}$ and $\\text{O}$ atoms)

            *Pass Criterion:* An orally active drug typically has **no more than 1 violation** of these criteria.
            - **Caffeine (`CN1C=NC2=C1C(=O)N(C(=O)N2C)C`)**: MW = 194.19, LogP = -1.03, HBD = 0, HBA = 6 → **0 Violations (Pass)**.
            """
        )


# -----------------------------------------------------------------------------
# PAGE 3: MOLECULAR VISUALIZATION
# -----------------------------------------------------------------------------

def render_visualization() -> None:
    """Render the Molecular Visualization page with 2D | 3D tabs and example buttons."""
    st.header("🔬 Molecular Visualization (2D & 3D)")
    st.write(
        "Enter any SMILES string or click a benchmark molecule below to inspect its "
        "2D structural formula and interactive 3D force-field minimized conformation."
    )

    st.caption("Quick-Load Example Molecules:")
    btn_cols = st.columns(len(EXAMPLE_MOLECULES))
    for col, (name, smi) in zip(btn_cols, EXAMPLE_MOLECULES.items()):
        with col:
            if st.button(
                name,
                key=f"viz_btn_{name}",
                use_container_width=True,
                help=f"Load {name} ({smi}) into the visualizer.",
            ):
                st.session_state["viz_smiles"] = smi
                st.toast(f"Loaded {name} ({smi})", icon="🧪")

    smiles_input = st.text_input(
        "SMILES Notation Input",
        value=st.session_state["viz_smiles"],
        placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O",
        help="Type a valid SMILES string (e.g., CCO for Ethanol or c1ccccc1 for Benzene).",
    )
    st.session_state["viz_smiles"] = smiles_input

    mol, err = validate_smiles(smiles_input)
    if err:
        st.error(f"⚠️ **Invalid SMILES Input:** {err}")
        return

    props = calculate_molecular_properties(smiles_input)
    if props:
        record_studied_molecule(props)

    tab_2d, tab_3d = st.tabs(["🖼️ 2D Depiction (RDKit)", "🌐 Interactive 3D Model (py3Dmol)"])

    with tab_2d:
        with st.spinner("Generating 2D depiction with RDKit..."):
            img_bytes = render_2d_image(smiles_input)
        if img_bytes:
            c1, c2 = st.columns([2, 1])
            with c1:
                st.image(img_bytes, caption=f"2D Structure — {props['Name']} ({props['Formula']})")
            with c2:
                st.markdown(
                    f"""
                    <div class="lab-card">
                        <h4>Structure Summary</h4>
                        <p><strong>Molecule:</strong> {props['Name']}</p>
                        <p><strong>Canonical SMILES:</strong> <code>{props['SMILES']}</code></p>
                        <p><strong>Formula:</strong> {props['Formula']}</p>
                        <p><strong>Molecular Weight:</strong> {props['Molecular Weight']} g/mol</p>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
            st.toast("2D structure rendered!", icon="✅")

    with tab_3d:
        ctrl1, ctrl2 = st.columns(2)
        with ctrl1:
            style_choice = st.selectbox(
                "3D Representation Style",
                options=["stick", "sphere", "line"],
                index=0,
                help="Choose 'stick' for covalent bonds, 'sphere' for CPK van der Waals radii, or 'line' for wireframe.",
            )
        with ctrl2:
            color_scheme = st.selectbox(
                "Atom Color Scheme",
                options=["default (CPK)", "cyanCarbon", "greenCarbon", "magentaCarbon", "yellowCarbon"],
                index=0,
                help="Select the carbon atom coloring scheme while keeping N (blue), O (red), and H (white).",
            )

        with st.spinner("Adding explicit hydrogens & optimizing 3D geometry (MMFF94)..."):
            molblock, err_3d = build_3d_molblock(smiles_input)

        if err_3d or not molblock:
            st.error(f"⚠️ {err_3d}")
        else:
            scheme_arg = "default" if color_scheme.startswith("default") else color_scheme
            view = py3Dmol.view(width=700, height=440)
            view.addModel(molblock, "mol")
            view.setStyle({style_choice: {"colorscheme": scheme_arg}})
            view.setBackgroundColor("#0F172A" if st.session_state["dark_mode"] else "#FFFFFF")
            view.zoomTo()
            render_3d_viewer(view, height=450, width=700)
            st.caption("🖱️ **Controls:** Left-click + drag to rotate • Scroll to zoom • Right-click + drag to pan.")


# -----------------------------------------------------------------------------
# PAGE 4: MOLECULAR PROPERTIES REPORT
# -----------------------------------------------------------------------------

def render_properties_report() -> None:
    """Render the Molecular Properties Report page with metrics, table, CSV & PDF exports."""
    st.header("📋 Molecular Properties Report")
    st.write(
        "Compute key physicochemical descriptors and check Lipinski's Rule of Five compliance "
        "for any organic molecule."
    )

    col_inp, col_preset = st.columns([2, 1])
    with col_preset:
        preset = st.selectbox(
            "Or Load an Example Molecule",
            options=["Custom Input"] + list(EXAMPLE_MOLECULES.keys()),
            help="Quickly populate the SMILES box with one of the reference molecules.",
        )
        if preset != "Custom Input":
            st.session_state["prop_smiles"] = EXAMPLE_MOLECULES[preset]

    with col_inp:
        smiles_input = st.text_input(
            "Enter Molecule SMILES",
            value=st.session_state["prop_smiles"],
            placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O",
            help="Enter a valid SMILES notation to calculate 9 physicochemical and Lipinski descriptors.",
        )
        st.session_state["prop_smiles"] = smiles_input

    _, err = validate_smiles(smiles_input)
    if err:
        st.error(f"⚠️ **Cannot Compute Properties:** {err}")
        return

    with st.spinner("Calculating molecular descriptors..."):
        props = calculate_molecular_properties(smiles_input)

    if not props:
        st.error("⚠️ Could not compute descriptors for the given input.")
        return

    record_studied_molecule(props)
    st.toast(f"Computed descriptors for {props['Name']}", icon="📊")

    badge_class = "badge-pass" if props["Lipinski Check"] == "Pass" else "badge-fail"
    st.markdown(
        f"**Molecule:** `{props['Name']}` &nbsp;|&nbsp; **Lipinski Rule of Five Status:** "
        f"<span class='{badge_class}'>Lipinski {props['Lipinski Check']} "
        f"({props['Lipinski Violations']} Violations)</span>",
        unsafe_allow_html=True,
    )
    st.write("")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Molecular Weight", f"{props['Molecular Weight']} g/mol", help="Target <= 500 Da for Lipinski Ro5")
    m2.metric("Chemical Formula", props["Formula"], help="Hill system molecular formula")
    m3.metric("LogP (Lipophilicity)", f"{props['LogP']}", help="Target <= 5.0 for Lipinski Ro5")
    m4.metric("TPSA", f"{props['TPSA']} Å²", help="Topological Polar Surface Area (typically < 140 Å² for oral drugs)")

    m5, m6, m7, m8 = st.columns(4)
    m5.metric("H-Bond Donors (HBD)", props["H-bond Donors"], help="Target <= 5 OH/NH groups")
    m6.metric("H-Bond Acceptors (HBA)", props["H-bond Acceptors"], help="Target <= 10 N/O atoms")
    m7.metric("Rotatable Bonds", props["Rotatable Bonds"], help="Measures conformational flexibility")
    m8.metric("Ring Count", props["Ring Count"], help="Total number of rings in the smallest set of smallest rings")

    st.subheader("Detailed Descriptor Table")
    df_single = pd.DataFrame([props])
    st.dataframe(df_single, use_container_width=True, hide_index=True)

    dl1, dl2 = st.columns(2)
    with dl1:
        csv_bytes = df_single.to_csv(index=False).encode("utf-8")
        st.download_button(
            "📥 Download Properties (CSV)",
            data=csv_bytes,
            file_name=f"molecular_properties_{props['Name']}.csv",
            mime="text/csv",
            use_container_width=True,
            help="Download this molecule's calculated descriptors as a comma-separated values (CSV) spreadsheet.",
        )
    with dl2:
        pdf_bytes = generate_pdf_report(
            title="Molecular Properties Report",
            student_name=st.session_state["student_name"],
            reg_number=st.session_state["reg_number"],
            properties_list=[props],
        )
        st.download_button(
            "📄 Download Properties Report (PDF)",
            data=pdf_bytes,
            file_name=f"molecular_properties_{props['Name']}.pdf",
            mime="application/pdf",
            use_container_width=True,
            help="Download a formatted academic PDF report containing student info and molecular descriptors.",
        )


# -----------------------------------------------------------------------------
# PAGE 5: STRUCTURE-PROPERTY ANALYSIS
# -----------------------------------------------------------------------------

def render_structure_property_analysis() -> None:
    """Render the Structure-Property Analysis page for multi-SMILES comparison."""
    st.header("📈 Structure-Property Relationship (SAR) Analysis")
    st.write(
        "Enter multiple SMILES strings (one per line) to compare physicochemical descriptors "
        "across a chemical series and visualize trends."
    )

    raw_lines = st.text_area(
        "Batch SMILES Input (One SMILES string per line)",
        value=st.session_state["batch_smiles"],
        height=155,
        help="Enter at least 2 SMILES strings separated by newlines to compare their molecular properties.",
    )
    st.session_state["batch_smiles"] = raw_lines

    lines = [line.strip() for line in raw_lines.splitlines() if line.strip()]
    if not lines:
        st.warning("⚠️ Please enter at least one SMILES string above. Suggested fix: Paste `CCO` and `c1ccccc1` on separate lines.")
        return

    valid_props: List[Dict[str, Any]] = []
    invalid_lines: List[str] = []

    with st.spinner("Analyzing molecular series..."):
        for smi in lines:
            props = calculate_molecular_properties(smi)
            if props:
                valid_props.append(props)
                record_studied_molecule(props)
            else:
                invalid_lines.append(smi)

    if invalid_lines:
        st.warning(
            f"⚠️ Skipped {len(invalid_lines)} invalid SMILES line(s): `{', '.join(invalid_lines)}`. "
            "Suggested fix: Verify ring closure numbers and atomic capitalization."
        )

    if not valid_props:
        st.error("⚠️ No valid SMILES strings found to plot.")
        return

    df_comp = pd.DataFrame(valid_props)
    st.subheader("Multi-Molecule Comparison Table")
    st.dataframe(df_comp, use_container_width=True, hide_index=True)

    col_bar, col_scatter = st.columns(2)

    with col_bar:
        metric_col = st.selectbox(
            "Select Descriptor for Bar Chart Comparison",
            options=[
                "Molecular Weight",
                "LogP",
                "TPSA",
                "H-bond Donors",
                "H-bond Acceptors",
                "Rotatable Bonds",
                "Ring Count",
            ],
            index=0,
            help="Choose which molecular descriptor to compare across all valid compounds.",
        )
        fig_bar = px.bar(
            df_comp,
            x="Name",
            y=metric_col,
            color="Lipinski Check",
            color_discrete_map={"Pass": "#14B8A6", "Fail": "#EF4444"},
            text_auto=".2f",
            title=f"{metric_col} Across Studied Molecules",
        )
        fig_bar.update_layout(xaxis_title="Molecule", yaxis_title=metric_col)
        st.plotly_chart(fig_bar, use_container_width=True)

    with col_scatter:
        st.write("")
        st.write("")
        fig_scatter = px.scatter(
            df_comp,
            x="Molecular Weight",
            y="LogP",
            color="Lipinski Check",
            size="TPSA",
            hover_name="Name",
            hover_data=["SMILES", "Formula", "H-bond Donors", "H-bond Acceptors"],
            color_discrete_map={"Pass": "#1E3A8A", "Fail": "#EF4444"},
            title="Chemical Space: Molecular Weight vs. LogP (Bubble Size = TPSA)",
        )
        fig_scatter.add_hline(y=5.0, line_dash="dash", line_color="#EF4444", annotation_text="Max LogP (5.0)")
        fig_scatter.add_vline(x=500.0, line_dash="dash", line_color="#EF4444", annotation_text="Max MW (500 Da)")
        st.plotly_chart(fig_scatter, use_container_width=True)


# -----------------------------------------------------------------------------
# PAGE 6: ASSESSMENT
# -----------------------------------------------------------------------------

def render_assessment() -> None:
    """Render 10-question multiple-choice quiz one question at a time with instant feedback."""
    st.header("📝 Cheminformatics Self-Assessment")
    st.write(
        "Answer all 10 multiple-choice questions to evaluate your mastery of SMILES, "
        "molecular graphs, descriptors, and Lipinski's Rule of Five."
    )

    total_q = len(QUIZ_QUESTIONS)
    q_idx = st.session_state["quiz_index"]
    answered_count = len(st.session_state["quiz_submitted"])

    st.progress(
        answered_count / total_q,
        text=f"Quiz Progress: {answered_count} / {total_q} Questions Answered",
    )

    if answered_count == total_q and q_idx >= total_q:
        score = sum(
            1
            for i, q in enumerate(QUIZ_QUESTIONS)
            if st.session_state["quiz_answers"].get(i) == q["answer"]
        )
        pct = score / total_q * 100
        passed = pct >= 60
        badge_cls = "badge-pass" if passed else "badge-fail"
        badge_txt = "🏆 Certified Cheminformatics Explorer (PASS)" if passed else "📚 Review Recommended (Keep Practicing)"

        st.markdown(
            f"""
            <div class="lab-card" style="text-align:center; padding: 2rem;">
                <h3>Final Assessment Score: {score} / {total_q} ({pct:.0f}%)</h3>
                <p style="margin: 1rem 0;"><span class="{badge_cls}">{badge_txt}</span></p>
                <p>Your score has been recorded in your session profile for the Final Project Report.</p>
            </div>
            """,
            unsafe_allow_html=True,
        )

        if st.button("🔄 Retake Assessment", help="Reset all quiz answers and start from Question 1."):
            st.session_state["quiz_index"] = 0
            st.session_state["quiz_answers"] = {}
            st.session_state["quiz_submitted"] = {}
            st.rerun()
        return

    q_idx = min(q_idx, total_q - 1)
    q_data = QUIZ_QUESTIONS[q_idx]

    st.subheader(f"Question {q_idx + 1} of {total_q}")
    st.markdown(f"**{q_data['q']}**")

    prev_choice = st.session_state["quiz_answers"].get(q_idx, q_data["options"][0])
    default_opt_idx = q_data["options"].index(prev_choice) if prev_choice in q_data["options"] else 0

    selected_opt = st.radio(
        "Select the best answer:",
        options=q_data["options"],
        index=default_opt_idx,
        key=f"radio_q_{q_idx}",
        disabled=st.session_state["quiz_submitted"].get(q_idx, False),
        help="Choose one option and click 'Check Answer' to receive instant feedback.",
    )

    c_sub, c_prev, c_next = st.columns([2, 1, 1])
    with c_sub:
        if not st.session_state["quiz_submitted"].get(q_idx, False):
            if st.button(
                "✔️ Check Answer",
                type="primary",
                use_container_width=True,
                help="Submit your choice for this question and view the explanation.",
            ):
                st.session_state["quiz_answers"][q_idx] = selected_opt
                st.session_state["quiz_submitted"][q_idx] = True
                st.toast("Answer recorded!", icon="📝")
                st.rerun()

    if st.session_state["quiz_submitted"].get(q_idx, False):
        user_ans = st.session_state["quiz_answers"].get(q_idx)
        if user_ans == q_data["answer"]:
            st.success(f"✅ **Correct!** {q_data['explanation']}")
        else:
            st.error(
                f"❌ **Incorrect.** You chose *{user_ans}*, but the correct answer is **{q_data['answer']}**.\n\n"
                f"**Explanation:** {q_data['explanation']}"
            )

    with c_prev:
        if st.button(
            "⬅️ Previous",
            disabled=(q_idx == 0),
            use_container_width=True,
            help="Navigate to the previous question.",
        ):
            st.session_state["quiz_index"] = max(0, q_idx - 1)
            st.rerun()

    with c_next:
        next_label = "🏁 View Final Score" if (q_idx == total_q - 1 and len(st.session_state["quiz_submitted"]) == total_q) else "Next ➡️"
        if st.button(
            next_label,
            disabled=(q_idx == total_q - 1 and len(st.session_state["quiz_submitted"]) < total_q),
            use_container_width=True,
            help="Advance to the next question or view your final score screen.",
        ):
            st.session_state["quiz_index"] = q_idx + 1
            st.rerun()


# -----------------------------------------------------------------------------
# PAGE 7: FINAL PROJECT REPORT
# -----------------------------------------------------------------------------

def render_final_report() -> None:
    """Render the Final Project Report summarizing student info, molecules, and quiz score."""
    st.header("🎓 Final Project Report")
    st.write(
        "Review your complete laboratory portfolio below and download the official PDF "
        "submission for your course instructor."
    )

    student_name = st.session_state["student_name"] or "Not Entered (Update in Sidebar)"
    reg_number = st.session_state["reg_number"] or "Not Entered (Update in Sidebar)"

    total_q = len(QUIZ_QUESTIONS)
    score = sum(
        1
        for i, q in enumerate(QUIZ_QUESTIONS)
        if st.session_state["quiz_answers"].get(i) == q["answer"]
    )
    answered_q = len(st.session_state["quiz_submitted"])

    c1, c2, c3 = st.columns(3)
    c1.metric("Student Name", student_name, help="Pulled automatically from the sidebar.")
    c2.metric("Registration Number", reg_number, help="Pulled automatically from the sidebar.")
    c3.metric(
        "Assessment Score",
        f"{score} / {total_q} ({score / total_q * 100:.0f}%)",
        delta=f"{answered_q}/{total_q} answered",
        help="Your current score from the Assessment module.",
    )

    if not st.session_state["studied_molecules"]:
        for smi in list(EXAMPLE_MOLECULES.values())[:3]:
            props = calculate_molecular_properties(smi)
            if props:
                record_studied_molecule(props)

    studied_list = list(st.session_state["studied_molecules"].values())

    st.subheader(f"Molecules Studied in This Session ({len(studied_list)})")
    df_studied = pd.DataFrame(studied_list)
    st.dataframe(df_studied, use_container_width=True, hide_index=True)

    st.subheader("Properties Summary Statistics")
    summary_cols = ["Molecular Weight", "LogP", "TPSA", "H-bond Donors", "H-bond Acceptors", "Rotatable Bonds"]
    st.dataframe(
        df_studied[summary_cols].describe().loc[["mean", "min", "max"]].round(2),
        use_container_width=True,
    )

    with st.spinner("Preparing downloadable PDF portfolio..."):
        pdf_bytes = generate_pdf_report(
            title="Final Project Portfolio",
            student_name=st.session_state["student_name"],
            reg_number=st.session_state["reg_number"],
            properties_list=studied_list,
            assessment_score=(score, total_q),
        )

    st.download_button(
        "📥 Download Final Project Report (PDF)",
        data=pdf_bytes,
        file_name=f"Cheminformatics_Final_Report_{st.session_state['reg_number'] or 'Student'}.pdf",
        mime="application/pdf",
        type="primary",
        use_container_width=True,
        help="Download your complete Cheminformatics Virtual Laboratory PDF report.",
    )


# -----------------------------------------------------------------------------
# MAIN APPLICATION ROUTER
# -----------------------------------------------------------------------------

def main() -> None:
    """Configure page, initialize state, render sidebar, and route to active page."""
    st.set_page_config(
        page_title="Cheminformatics Virtual Laboratory",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    init_session_state()
    active_page = render_sidebar()
    inject_custom_css(dark_mode=st.session_state["dark_mode"])

    page_router = {
        "Home": render_home,
        "Theory": render_theory,
        "Molecular Visualization": render_visualization,
        "Molecular Properties Report": render_properties_report,
        "Structure-Property Analysis": render_structure_property_analysis,
        "Assessment": render_assessment,
        "Final Project Report": render_final_report,
    }
    page_router[active_page]()


if __name__ == "__main__":
    main()

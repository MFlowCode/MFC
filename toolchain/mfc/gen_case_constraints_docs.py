#!/usr/bin/env python3
"""
Generate human-readable documentation for MFC case parameter constraints.

Parses toolchain/mfc/case_validator.py, extracts all `self.prohibit(...)` rules,
maps them to parameters and stages, and emits Markdown to stdout.

Also generates case design playbook from curated working examples.
"""  # pylint: disable=too-many-lines

import json
import sys
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Iterable, Any
from collections import defaultdict

HERE = Path(__file__).resolve().parent
CASE_VALIDATOR_PATH = HERE / "case_validator.py"
REPO_ROOT = HERE.parent.parent
EXAMPLES_DIR = REPO_ROOT / "examples"

# Make the params package importable
_toolchain_dir = str(HERE.parent)
if _toolchain_dir not in sys.path:
    sys.path.insert(0, _toolchain_dir)

from mfc.params import CONSTRAINTS, DEPENDENCIES, get_value_label  # noqa: E402  pylint: disable=wrong-import-position
from mfc.params.ast_analyzer import (  # noqa: E402  pylint: disable=wrong-import-position
    Rule, classify_message, feature_title,
    analyze_case_validator,
)


# ---------------------------------------------------------------------------
# Case Playbook Generation (from working examples)
# ---------------------------------------------------------------------------

@dataclass
class PlaybookEntry:
    """A curated example case for the playbook"""
    case_dir: str
    title: str
    description: str
    level: str  # "Beginner", "Intermediate", "Advanced"
    tags: List[str]


# Curated list of hero examples
PLAYBOOK_EXAMPLES = [
    PlaybookEntry(
        "2D_shockbubble",
        "2D Shock-Bubble Interaction",
        "Two-fluid shock-interface benchmark. Classic validation case for compressible multiphase flows.",
        "Beginner",
        ["2D", "Multiphase", "Shock"]
    ),
    PlaybookEntry(
        "1D_bubblescreen",
        "1D Bubble Screen",
        "Euler-Euler ensemble-averaged bubble dynamics through shock wave.",
        "Intermediate",
        ["1D", "Bubbles", "Euler-Euler"]
    ),
    PlaybookEntry(
        "2D_lagrange_bubblescreen",
        "2D Lagrangian Bubble Screen",
        "Individual bubble tracking with Euler-Lagrange method.",
        "Intermediate",
        ["2D", "Bubbles", "Euler-Lagrange"]
    ),
    PlaybookEntry(
        "2D_phasechange_bubble",
        "2D Phase Change Bubble",
        "Phase change and cavitation modeling with 6-equation model.",
        "Advanced",
        ["2D", "Phase-change", "Cavitation"]
    ),
    PlaybookEntry(
        "2D_orszag_tang",
        "2D Orszag-Tang MHD Vortex",
        "Magnetohydrodynamics test problem with complex vortex structures.",
        "Intermediate",
        ["2D", "MHD"]
    ),
    PlaybookEntry(
        "2D_ibm_airfoil",
        "2D IBM Airfoil",
        "Immersed boundary method around a NACA airfoil geometry.",
        "Intermediate",
        ["2D", "IBM", "Geometry"]
    ),
    PlaybookEntry(
        "2D_viscous_shock_tube",
        "2D Viscous Shock Tube",
        "Shock tube with viscous effects and heat transfer.",
        "Intermediate",
        ["2D", "Viscous", "Shock"]
    ),
    PlaybookEntry(
        "3D_TaylorGreenVortex",
        "3D Taylor-Green Vortex",
        "Classic 3D turbulence benchmark with viscous dissipation.",
        "Advanced",
        ["3D", "Viscous", "Turbulence"]
    ),
    PlaybookEntry(
        "2D_IGR_triple_point",
        "2D IGR Triple Point",
        "Triple point problem using Iterative Generalized Riemann solver.",
        "Advanced",
        ["2D", "IGR", "Multiphase"]
    ),
]


def validate_playbook_examples():
    """Check that all curated examples exist and error if any are missing"""
    missing = []
    for entry in PLAYBOOK_EXAMPLES:
        case_path = EXAMPLES_DIR / entry.case_dir / "case.py"
        if not case_path.exists():
            missing.append(entry.case_dir)

    if missing:
        print("=" * 70, file=sys.stderr)
        print("ERROR: Missing playbook examples:", file=sys.stderr)
        for example in missing:
            print(f"  - {example}", file=sys.stderr)
        print("\nPlease update PLAYBOOK_EXAMPLES in:", file=sys.stderr)
        print(f"  {Path(__file__).relative_to(REPO_ROOT)}", file=sys.stderr)
        print("\nRemove the missing examples from the list or restore them.", file=sys.stderr)
        print("=" * 70, file=sys.stderr)
        sys.exit(1)


def load_case_params(case_dir: str) -> Dict[str, Any]:
    """Load parameters from a case.py file"""
    case_path = EXAMPLES_DIR / case_dir / "case.py"
    if not case_path.exists():
        return {}

    try:
        result = subprocess.run(
            ["python3", str(case_path)],
            capture_output=True,
            text=True,
            timeout=10,
            check=True
        )
        params = json.loads(result.stdout)
        return params
    except (subprocess.CalledProcessError, json.JSONDecodeError, subprocess.TimeoutExpired) as e:
        print(f"WARNING: Failed to load params from {case_path}: {e}", file=sys.stderr)
        return {}


def summarize_case_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """Extract key features from case parameters"""
    return {
        "model_eqns": params.get("model_eqns"),
        "num_fluids": params.get("num_fluids"),
        "surface_tension": params.get("surface_tension") == "T",
        "bubbles_euler": params.get("bubbles_euler") == "T",
        "bubbles_lagrange": params.get("bubbles_lagrange") == "T",
        "qbmm": params.get("qbmm") == "T",
        "polydisperse": params.get("polydisperse") == "T",
        "mhd": params.get("mhd") == "T",
        "relax": params.get("relax") == "T",
        "hypoelasticity": params.get("hypoelasticity") == "T",
        "viscous": params.get("viscous") == "T",
        "ib": params.get("ib") == "T",
        "igr": params.get("igr") == "T",
        "acoustic_source": params.get("acoustic_source") == "T",
        "cyl_coord": params.get("cyl_coord") == "T",
        "m": params.get("m"),
        "n": params.get("n", 0),
        "p": params.get("p", 0),
        "recon_type": params.get("recon_type", 1),
        "weno_order": params.get("weno_order"),
        "muscl_order": params.get("muscl_order"),
        "riemann_solver": params.get("riemann_solver"),
        "time_stepper": params.get("time_stepper"),
    }


def get_model_name(model_eqns: int | None) -> str:
    """Get human-friendly model name from schema."""
    if model_eqns is None:
        return "Not specified"
    return get_value_label("model_eqns", model_eqns) or "Not specified"


def get_riemann_solver_name(solver: int | None) -> str:
    """Get Riemann solver name from schema."""
    if solver is None:
        return "Not specified"
    return get_value_label("riemann_solver", solver) or "Not specified"


def get_time_stepper_name(stepper: int | None) -> str:
    """Get time stepper name from schema."""
    if stepper is None:
        return "Not specified"
    return get_value_label("time_stepper", stepper) or "Not specified"


def render_playbook_card(entry: PlaybookEntry, summary: Dict[str, Any]) -> str:  # pylint: disable=too-many-branches,too-many-statements
    """Render a single playbook entry as Markdown"""
    lines = []

    tags_str = " · ".join(entry.tags)
    level_emoji = {"Beginner": "🟢", "Intermediate": "🟡", "Advanced": "🔴"}.get(entry.level, "")

    lines.append("<details>")
    lines.append(f'<summary><b>{entry.title}</b> {level_emoji} <i>{entry.level}</i> · <code>{entry.case_dir}</code></summary>\n')
    lines.append(f"**{entry.description}**\n")
    lines.append(f"**Tags:** {tags_str}\n")

    lines.append("**Physics Configuration:**\n")
    lines.append(f"- **Model:** {get_model_name(summary['model_eqns'])} (`model_eqns = {summary['model_eqns']}`)")

    if summary['num_fluids'] is not None:
        lines.append(f"- **Number of fluids:** {summary['num_fluids']}")

    # Dimensionality
    n, p = summary['n'], summary['p']
    dim_str = "3D" if p > 0 else ("2D" if n > 0 else "1D")
    lines.append(f"- **Dimensionality:** {dim_str}")

    if summary['cyl_coord']:
        lines.append("- **Coordinates:** Cylindrical/Axisymmetric")

    # Active features
    active_features = []
    if summary['bubbles_euler']:
        active_features.append("Euler-Euler bubbles")
    if summary['bubbles_lagrange']:
        active_features.append("Euler-Lagrange bubbles")
    if summary['qbmm']:
        active_features.append("QBMM")
    if summary['polydisperse']:
        active_features.append("Polydisperse")
    if summary['surface_tension']:
        active_features.append("Surface tension")
    if summary['mhd']:
        active_features.append("MHD")
    if summary['relax']:
        active_features.append("Phase change")
    if summary['hypoelasticity']:
        active_features.append("Hypoelasticity")
    if summary['viscous']:
        active_features.append("Viscous")
    if summary['ib']:
        active_features.append("Immersed boundaries")
    if summary['igr']:
        active_features.append("IGR solver")
    if summary['acoustic_source']:
        active_features.append("Acoustic sources")

    if active_features:
        lines.append(f"- **Active features:** {', '.join(active_features)}")

    # Numerics
    lines.append("\n**Numerical Methods:**\n")

    if summary['recon_type'] == 1 and summary['weno_order']:
        lines.append(f"- **Reconstruction:** WENO-{summary['weno_order']}")
    elif summary['recon_type'] == 2 and summary['muscl_order']:
        lines.append(f"- **Reconstruction:** MUSCL (order {summary['muscl_order']})")

    if summary['riemann_solver']:
        solver_name = get_riemann_solver_name(summary['riemann_solver'])
        lines.append(f"- **Riemann solver:** {solver_name} (`riemann_solver = {summary['riemann_solver']}`)")

    if summary['time_stepper']:
        stepper_name = get_time_stepper_name(summary['time_stepper'])
        lines.append(f"- **Time stepping:** {stepper_name}")

    # Links
    lines.append("\n**Related Documentation:**")
    lines.append(f"- [Model Equations (model_eqns = {summary['model_eqns']})](#model-equations)")

    if summary['riemann_solver']:
        lines.append("- [Riemann Solvers](#riemann-solvers)")

    if summary['bubbles_euler'] or summary['bubbles_lagrange']:
        lines.append("- [Bubble Models](#bubble-models)")

    if summary['mhd']:
        lines.append("- [MHD](#compat-physics-models)")

    if summary['ib']:
        lines.append("- [Immersed Boundaries](#compat-geometry)")

    if summary['viscous']:
        lines.append("- [Viscosity](#compat-physics-models)")

    lines.append("\n</details>\n")
    return "\n".join(lines)


def generate_playbook() -> str:
    """Generate complete playbook from curated examples"""
    lines = []

    # Validate examples - will exit(1) if any are missing
    validate_playbook_examples()

    lines.append("## 🧩 Case Design Playbook {#case-design-playbook}\n")
    lines.append(
        "> **Learn by example:** The cases below are curated from MFC's `examples/` "
        "directory and are validated, working configurations. "
        "Use them as blueprints for building your own simulations.\n"
    )

    # Group by level
    for level in ["Beginner", "Intermediate", "Advanced"]:
        level_entries = [e for e in PLAYBOOK_EXAMPLES if e.level == level]
        if not level_entries:
            continue

        level_emoji = {"Beginner": "🟢", "Intermediate": "🟡", "Advanced": "🔴"}.get(level, "")
        lines.append(f"\n### {level_emoji} {level} Examples\n")

        for entry in level_entries:
            try:
                params = load_case_params(entry.case_dir)
                if not params:
                    continue
                summary = summarize_case_params(params)
                card = render_playbook_card(entry, summary)
                lines.append(card)
            except Exception as e:  # pylint: disable=broad-except
                print(f"WARNING: Failed to process playbook entry '{entry.case_dir}': {e}", file=sys.stderr)
                continue

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(rules: Iterable[Rule]) -> str:  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """
    Render user-friendly compatibility tables and summaries.
    """
    # Group by parameter
    by_param: Dict[str, List[Rule]] = defaultdict(list)
    for r in rules:
        if not r.params:
            continue
        for p in r.params:
            by_param[p].append(r)

    lines: List[str] = []

    lines.append("@page case_constraints Case Creator Guide\n")
    lines.append("# Case Creator Guide\n")
    lines.append(
        "> **Quick reference** for building MFC cases: working examples, compatibility rules, "
        "and configuration requirements.\n"
    )
    lines.append(
        "> Auto-generated from `case_validator.py` and `examples/`.\n"
    )

    # Add playbook at the top
    playbook = generate_playbook()
    lines.append(playbook)

    # Define major feature groups (excluding IGR)
    major_features = {
        "Physics Models": ["mhd", "surface_tension", "hypoelasticity", "hyperelasticity", "relax", "viscous", "acoustic_source"],
        "Bubble Models": ["bubbles_euler", "bubbles_lagrange", "qbmm", "polydisperse", "adv_n"],
        "Numerics": ["riemann_solver", "weno_order", "muscl_order"],
        "Geometry": ["ib", "cyl_coord"],
    }

    # 1. Quick Start: Common Configurations
    lines.append("## 🚀 Common Configuration Patterns {#common-configuration-patterns}\n")
    lines.append("Start with these proven combinations:\n")
    lines.append("")
    lines.append("<details open>")
    lines.append("<summary><b>💧 Multiphase Flow (Bubbles)</b></summary>\n")
    lines.append("```python")
    lines.append("'model_eqns': 2,              # 5-equation model")
    lines.append("'num_fluids': 2,              # Two-phase")
    lines.append("'bubbles_euler': 'T',         # Ensemble-averaged bubbles")
    lines.append("'riemann_solver': 2,          # HLLC")
    lines.append("'avg_state': 2,               # Arithmetic average")
    lines.append("```")
    lines.append("</details>\n")

    lines.append("<details>")
    lines.append("<summary><b>⚡ Magnetohydrodynamics (MHD)</b></summary>\n")
    lines.append("```python")
    lines.append("'model_eqns': 2,              # 5-equation model")
    lines.append("'num_fluids': 1,              # Single component")
    lines.append("'mhd': 'T',                   # Enable MHD")
    lines.append("'riemann_solver': 1,          # HLL (or 4 for HLLD)")
    lines.append("```")
    lines.append("</details>\n")

    lines.append("<details>")
    lines.append("<summary><b>🌡️ Phase Change</b></summary>\n")
    lines.append("```python")
    lines.append("'model_eqns': 3,              # 6-equation model")
    lines.append("'num_fluids': 2,              # Two-phase")
    lines.append("'relax': 'T',                 # Phase change")
    lines.append("'riemann_solver': 2,          # HLLC")
    lines.append("```")
    lines.append("</details>\n")

    lines.append("<details>")
    lines.append("<summary><b>💎 Elastic Materials</b></summary>\n")
    lines.append("```python")
    lines.append("'model_eqns': 2,              # 5-equation model")
    lines.append("'hypoelasticity': 'T',        # Elastic solids")
    lines.append("'riemann_solver': 1,          # HLL")
    lines.append("```")
    lines.append("</details>\n")

    # 2. Feature Compatibility Matrix (simplified, no IGR column)
    lines.append("## 📊 Feature Compatibility {#feature-compatibility}\n")
    lines.append("What works together:\n")

    for category, features in major_features.items():  # pylint: disable=too-many-nested-blocks
        cat_id = "compat-" + category.lower().replace(" ", "-")
        lines.append(f"\n### {category} {{#{cat_id}}}\n")

        # Build compatibility info (exclude IGR from incompatibilities)
        compat_info = {}
        for feat in features:
            if feat not in by_param:
                continue

            incomp = []
            req = []
            for rule in by_param[feat]:
                kind = classify_message(rule.message)
                msg = rule.message

                # Skip IGR-related incompatibilities
                if "IGR" in msg:
                    continue

                if kind == "incompatibility":
                    for other_feat in features:
                        if other_feat != feat and other_feat in rule.params:
                            incomp.append(feature_title(other_feat))
                elif kind == "requirement":
                    for other_feat in features:
                        if other_feat != feat and other_feat in rule.params:
                            req.append(feature_title(other_feat))

            compat_info[feat] = {"incomp": list(set(incomp)), "req": list(set(req))}

        # Render as compact table
        lines.append("| Feature | Requirements | Notes |")
        lines.append("|---------|-------------|-------|")

        for feat in features:
            if feat not in by_param:
                continue
            title = feature_title(feat)
            info = compat_info.get(feat, {"incomp": [], "req": []})

            # Get key requirement
            reqs = []
            for rule in by_param[feat]:
                if classify_message(rule.message) == "requirement":
                    msg = rule.message
                    if "model_eqns" in msg:
                        reqs.append("Specific model")
                        break

            req_str = reqs[0] if reqs else "—"

            # Get short note
            incomp_with = [i for i in info["incomp"] if i not in ["IGR"]][:2]
            if incomp_with:
                note = f"⚠️ Not with: {', '.join(incomp_with)}"
            else:
                note = "✓ General use"

            lines.append(f"| {title} | {req_str} | {note} |")

        lines.append("")

    # 3. Model Equations (data-driven from schema)
    lines.append("## 🔢 Model Equations {#model-equations}\n")
    lines.append("Choose your governing equations:\n")
    lines.append("")

    def _format_model_requirements(val: int) -> str:
        """Auto-generate requirements string from DEPENDENCIES['model_eqns']['when_value']."""
        me_dep = DEPENDENCIES.get("model_eqns", {})
        wv = me_dep.get("when_value", {}).get(val, {})
        if not wv:
            return ""
        parts = []
        if "requires" in wv:
            parts.extend(f"Set `{r}`" for r in wv["requires"])
        if "requires_value" in wv:
            for rv_param, rv_vals in wv["requires_value"].items():
                labeled = [f"`{v}` ({get_value_label(rv_param, v)})" for v in rv_vals]
                parts.append(f"`{rv_param}` = {' or '.join(labeled)}")
        return ", ".join(parts)

    # Curated editorial notes keyed by model_eqns value
    _model_notes = {
        1: {
            "use_for": "Single-fluid compressible flow",
            "note": "Cannot use `num_fluids`, bubbles, or certain WENO variants",
        },
        2: {
            "use_for": "Multiphase, bubbles, elastic materials, MHD",
            "note": "Compatible with most physics models",
        },
        3: {
            "use_for": "Phase change, cavitation",
            "note": "Not compatible with bubbles or 3D cylindrical",
        },
        4: {
            "use_for": "Single-component flows with bubbles",
        },
    }

    # Auto-populate requirements from schema
    for _val, _note in _model_notes.items():
        _auto_reqs = _format_model_requirements(_val)
        if _auto_reqs:
            _note["requirements"] = _auto_reqs

    model_constraint = CONSTRAINTS["model_eqns"]
    for val in model_constraint["choices"]:
        label = get_value_label("model_eqns", val)
        notes = _model_notes.get(val, {})
        lines.append("<details>")
        lines.append(f"<summary><b>Model {val}: {label}</b></summary>\n")
        if "use_for" in notes:
            lines.append(f"- **Use for:** {notes['use_for']}")
        lines.append(f"- **Value:** `model_eqns = {val}`")
        if "requirements" in notes:
            lines.append(f"- **Requirements:** {notes['requirements']}")
        if "note" in notes:
            lines.append(f"- **Note:** {notes['note']}")
        lines.append("</details>\n")

    # 4. Riemann Solvers (data-driven from schema)
    # Curated editorial notes keyed by riemann_solver value
    _solver_notes = {
        1: {"best_for": "MHD, elastic materials", "requirements": "—"},
        2: {"best_for": "Bubbles, phase change, multiphase", "requirements": "`avg_state=2` for bubbles"},
        3: {"best_for": "High accuracy (expensive)", "requirements": "—"},
        4: {"best_for": "MHD (advanced)", "requirements": "MHD only, no relativity"},
        5: {"best_for": "Robust fallback", "requirements": "Not with cylindrical+viscous"},
    }

    lines.append("## ⚙️ Riemann Solvers {#riemann-solvers}\n")
    lines.append("| Solver | `riemann_solver` | Best For | Requirements |")
    lines.append("|--------|-----------------|----------|-------------|")

    solver_constraint = CONSTRAINTS["riemann_solver"]
    for val in solver_constraint["choices"]:
        label = get_value_label("riemann_solver", val)
        notes = _solver_notes.get(val, {})
        best = notes.get("best_for", "—")
        reqs = notes.get("requirements", "—")
        lines.append(f"| **{label}** | `{val}` | {best} | {reqs} |")

    lines.append("")

    # 5. Bubble Models (data-driven from schema dependencies + curated notes)
    if "bubbles_euler" in by_param or "bubbles_lagrange" in by_param:
        lines.append("## 💧 Bubble Models {#bubble-models}\n")
        lines.append("")

        # Euler-Euler: inject schema dependency info (data-driven)
        lines.append("<details>")
        lines.append("<summary><b>Euler-Euler (`bubbles_euler`)</b></summary>\n")
        lines.append("**Requirements:**")
        be_dep = DEPENDENCIES.get("bubbles_euler", {})
        be_when_true = be_dep.get("when_true", {})
        be_rv = be_when_true.get("requires_value", {})
        for rv_param, rv_vals in be_rv.items():
            labeled = [f"`{v}` ({get_value_label(rv_param, v)})" for v in rv_vals]
            lines.append(f"- `{rv_param}` = {' or '.join(labeled)}")
        be_recs = be_when_true.get("recommends", [])
        if be_recs:
            lines.append(f"- Recommended to also set: {', '.join(f'`{r}`' for r in be_recs)}")
        lines.append("")
        lines.append("**Extensions:**")
        # Inject polydisperse dependency
        pd_dep = DEPENDENCIES.get("polydisperse", {})
        pd_reqs = pd_dep.get("when_true", {}).get("requires", [])
        pd_req_str = f" (requires {', '.join(f'`{r}`' for r in pd_reqs)})" if pd_reqs else ""
        lines.append(f"- `polydisperse = T`: Multiple bubble sizes{pd_req_str}, odd `nb > 1`")
        # Inject qbmm dependency
        qb_dep = DEPENDENCIES.get("qbmm", {})
        qb_recs = qb_dep.get("when_true", {}).get("recommends", [])
        qb_rec_str = f" (recommends {', '.join(f'`{r}`' for r in qb_recs)})" if qb_recs else ""
        lines.append(f"- `qbmm = T`: Quadrature method{qb_rec_str}, requires `nnode = 4`")
        lines.append("- `adv_n = T`: Number density advection (requires `num_fluids = 1`)")
        lines.append("</details>\n")

        lines.append("<details>")
        lines.append("<summary><b>Euler-Lagrange (`bubbles_lagrange`)</b></summary>\n")
        lines.append("**Requirements:**")
        lines.append("- `n > 0` (2D or 3D only)")
        lines.append("- `file_per_process = F`")
        lines.append("- Not compatible with `model_eqns = 3`\n")
        lines.append("**Note:** Tracks individual bubbles")
        lines.append("</details>\n")

    # 6. Condensed Parameter Reference (auto-collected from schema)
    lines.append("## 📖 Quick Parameter Reference {#quick-parameter-reference}\n")
    lines.append("Key parameters and their constraints:\n")

    # Auto-collect all params that have CONSTRAINTS or DEPENDENCIES entries
    quick_ref_params = sorted(set(CONSTRAINTS.keys()) | set(DEPENDENCIES.keys()))

    for param in quick_ref_params:
        title = feature_title(param)

        # Gather schema info
        constraint = CONSTRAINTS.get(param, {})
        dep = DEPENDENCIES.get(param, {})

        # Gather AST-extracted rules
        rules_for_param = by_param.get(param, [])
        requirements = []
        incompatibilities = []
        ranges = []
        warnings = []

        for rule in rules_for_param:
            msg = rule.message
            if "IGR" in msg:
                continue
            if rule.severity == "warning":
                warnings.append(msg)
                continue
            kind = classify_message(msg)
            if kind == "requirement":
                requirements.append(msg)
            elif kind == "incompatibility":
                incompatibilities.append(msg)
            elif kind == "range":
                ranges.append(msg)

        # Build schema constraint summary
        schema_parts = []
        if "choices" in constraint:
            labels = constraint.get("value_labels", {})
            if labels:
                items = [f"`{v}` = {labels[v]}" for v in constraint["choices"] if v in labels]
                schema_parts.append("Choices: " + ", ".join(items))
            else:
                schema_parts.append(f"Choices: {constraint['choices']}")
        if "min" in constraint:
            schema_parts.append(f"Min: {constraint['min']}")
        if "max" in constraint:
            schema_parts.append(f"Max: {constraint['max']}")

        # Build dependency summary
        dep_parts = []

        def _render_cond_parts(trigger_str, cond_dict):
            """Render a condition dict into dep_parts entries."""
            if "requires" in cond_dict:
                dep_parts.append(f"When {trigger_str}, requires: {', '.join(f'`{r}`' for r in cond_dict['requires'])}")
            if "requires_value" in cond_dict:
                rv_items = []
                for rv_p, rv_vs in cond_dict["requires_value"].items():
                    labeled = [f"`{v}` ({get_value_label(rv_p, v)})" for v in rv_vs]
                    rv_items.append(f"`{rv_p}` = {' or '.join(labeled)}")
                dep_parts.append(f"When {trigger_str}, requires {', '.join(rv_items)}")
            if "recommends" in cond_dict:
                dep_parts.append(f"When {trigger_str}, recommends: {', '.join(f'`{r}`' for r in cond_dict['recommends'])}")

        for cond_key in ["when_true", "when_set"]:
            cond = dep.get(cond_key, {})
            if cond:
                trigger = "enabled" if cond_key == "when_true" else "set"
                _render_cond_parts(trigger, cond)

        if "when_value" in dep:
            for wv_val, wv_cond in dep["when_value"].items():
                _render_cond_parts(f"= {wv_val}", wv_cond)

        # Skip if nothing to show
        if not (schema_parts or dep_parts or requirements or incompatibilities or ranges or warnings):
            continue

        lines.append(f"\n<details>")
        lines.append(f"<summary><b>{title}</b> (`{param}`)</summary>\n")

        if schema_parts:
            lines.append("**Schema constraints:**")
            for sp in schema_parts:
                lines.append(f"- {sp}")
            lines.append("")

        if dep_parts:
            lines.append("**Dependencies:**")
            for dp in dep_parts:
                lines.append(f"- {dp}")
            lines.append("")

        if requirements:
            lines.append("**Requirements** (errors):")
            for req in requirements[:3]:
                lines.append(f"- {req}")
            lines.append("")

        if incompatibilities:
            lines.append("**Incompatibilities** (errors):")
            for inc in incompatibilities[:3]:
                lines.append(f"- {inc}")
            lines.append("")

        if ranges:
            lines.append("**Valid values** (errors):")
            for rng in ranges[:2]:
                lines.append(f"- {rng}")
            lines.append("")

        if warnings:
            lines.append("**Warnings** (non-fatal):")
            for w in warnings[:3]:
                lines.append(f"- ⚠️ {w}")
            lines.append("")

        lines.append("</details>\n")

    # 7. Physics Warnings section (non-fatal checks from self.warn() calls)
    all_warnings = [r for r in rules if r.severity == "warning"]
    if all_warnings:
        lines.append("## ⚠️ Physics Warnings {#physics-warnings}\n")
        lines.append(
            "These checks are **non-fatal** — they print a yellow warning but do not abort the run. "
            "They catch common mistakes in initial conditions and EOS parameters.\n"
        )

        # Group by method
        warnings_by_method: Dict[str, List[Rule]] = defaultdict(list)
        for r in all_warnings:
            warnings_by_method[r.method].append(r)

        # Method name -> human-readable title
        method_titles = {
            "check_volume_fraction_sum": "Volume Fraction Sum",
            "check_alpha_rho_consistency": "Alpha-Rho Consistency",
            "check_eos_parameter_sanity": "EOS Parameter Sanity",
        }

        lines.append("| Check | Stage | Description |")
        lines.append("|-------|-------|-------------|")
        for method, method_rules in sorted(warnings_by_method.items()):
            title = method_titles.get(method, method.replace("check_", "").replace("_", " ").title())
            stages_str = ", ".join(sorted(method_rules[0].stages)) if method_rules[0].stages else "all"
            # Deduplicate messages (loop-expanded may repeat patterns)
            seen_msgs = set()
            descs = []
            for r in method_rules:
                if r.message not in seen_msgs:
                    seen_msgs.add(r.message)
                    descs.append(r.message)
            desc_str = "; ".join(descs[:2])
            if len(descs) > 2:
                desc_str += f" (+{len(descs)-2} more)"
            lines.append(f"| **{title}** | {stages_str} | {desc_str} |")

        lines.append("")

    # Add a footer with link to full validator
    lines.append("\n---\n")
    lines.append("💡 **Tip:** If you encounter a validation error, check the relevant section above or ")
    lines.append("review [`case_validator.py`](https://github.com/MFlowCode/MFC/blob/master/toolchain/mfc/case_validator.py) for complete validation logic.\n")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(as_string: bool = False) -> str:
    """Generate case constraints documentation. Returns markdown string."""
    analysis = analyze_case_validator(CASE_VALIDATOR_PATH)
    md = render_markdown(analysis["rules"])
    if not as_string:
        print(md)
    return md


if __name__ == "__main__":
    main()

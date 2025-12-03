#!/usr/bin/env python3
"""
Generate human-readable documentation for MFC case parameter constraints.

Parses toolchain/mfc/case_validator.py, extracts all `self.prohibit(...)` rules,
maps them to parameters and stages, and emits Markdown to stdout.

Also generates case design playbook from curated working examples.
"""

from __future__ import annotations

import ast
import json
import sys
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Set, Iterable, Any
from collections import defaultdict

HERE = Path(__file__).resolve().parent
CASE_VALIDATOR_PATH = HERE / "case_validator.py"
REPO_ROOT = HERE.parent.parent
EXAMPLES_DIR = REPO_ROOT / "examples"


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Rule:
    method: str                      # e.g. "check_igr_simulation"
    lineno: int                      # line number of the prohibit call
    params: List[str]                # case parameter names used in condition
    message: str                     # user-friendly error message
    stages: Set[str] = field(default_factory=set)  # e.g. {"simulation", "pre_process"}


# ---------------------------------------------------------------------------
# AST analysis: methods, call graph, rules
# ---------------------------------------------------------------------------

class CaseValidatorAnalyzer(ast.NodeVisitor):
    """
    Analyzes the CaseValidator class:

    - collects all methods
    - builds a call graph between methods
    - extracts all self.prohibit(...) rules
    """

    def __init__(self):
        super().__init__()
        self.in_case_validator = False
        self.current_method: str | None = None

        self.methods: Dict[str, ast.FunctionDef] = {}
        self.call_graph: Dict[str, Set[str]] = defaultdict(set)
        self.rules: List[Rule] = []

        # Stack of {local_name -> param_name} maps, one per method
        self.local_param_stack: List[Dict[str, str]] = []

    # --- top-level entrypoint ---

    def visit_ClassDef(self, node: ast.ClassDef):
        if node.name == "CaseValidator":
            self.in_case_validator = True
            # collect methods
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    self.methods[item.name] = item
            # now analyze all methods
            for method in self.methods.values():
                self._analyze_method(method)
            self.in_case_validator = False
        else:
            self.generic_visit(node)

    # --- per-method analysis ---

    def _analyze_method(self, func: ast.FunctionDef):
        """Analyze a single method: local param mapping, call graph, rules."""
        self.current_method = func.name
        local_param_map = self._build_local_param_map(func)
        self.local_param_stack.append(local_param_map)
        self.generic_visit(func)
        self.local_param_stack.pop()
        self.current_method = None

    def _build_local_param_map(self, func: ast.FunctionDef) -> Dict[str, str]:  # pylint: disable=too-many-nested-blocks
        """
        Look for assignments like:
            igr = self.get('igr', 'F') == 'T'
            model_eqns = self.get('model_eqns')
        and record local_name -> 'param_name'.
        """
        m: Dict[str, str] = {}
        for stmt in func.body:  # pylint: disable=too-many-nested-blocks
            if isinstance(stmt, ast.Assign):
                # Handle both direct calls and comparisons
                value = stmt.value
                # Unwrap comparisons like "self.get('igr', 'F') == 'T'"
                if isinstance(value, ast.Compare):
                    value = value.left

                if isinstance(value, ast.Call):
                    call = value
                    if (  # pylint: disable=too-many-boolean-expressions
                        isinstance(call.func, ast.Attribute)
                        and isinstance(call.func.value, ast.Name)
                        and call.func.value.id == "self"
                        and call.func.attr == "get"
                        and call.args
                        and isinstance(call.args[0], ast.Constant)
                        and isinstance(call.args[0].value, str)
                    ):
                        param_name = call.args[0].value
                        for target in stmt.targets:
                            if isinstance(target, ast.Name):
                                m[target.id] = param_name
        return m

    # --- visit calls to build call graph + rules ---

    def visit_Call(self, node: ast.Call):
        # record method call edges: self.some_method(...)
        if (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "self"
            and isinstance(node.func.attr, str)
        ):
            callee = node.func.attr
            if self.current_method is not None:
                # method call on self
                self.call_graph[self.current_method].add(callee)

        # detect self.prohibit(<condition>, "<message>")
        if (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "self"
            and node.func.attr == "prohibit"
            and len(node.args) >= 2
        ):
            condition, msg_node = node.args[0], node.args[1]
            if isinstance(msg_node, ast.Constant) and isinstance(msg_node.value, str):
                params = sorted(self._extract_params(condition))
                rule = Rule(
                    method=self.current_method or "<unknown>",
                    lineno=node.lineno,
                    params=params,
                    message=msg_node.value,
                )
                self.rules.append(rule)

        self.generic_visit(node)

    def _extract_params(self, condition: ast.AST) -> Set[str]:
        """
        Collect parameter names used in the condition via:
          - local variables mapped from self.get(...)
          - direct self.get('param_name', ...) calls
        """
        params: Set[str] = set()
        local_map = self.local_param_stack[-1] if self.local_param_stack else {}

        for node in ast.walk(condition):
            # local names
            if isinstance(node, ast.Name) and node.id in local_map:
                params.add(local_map[node.id])

            # direct self.get('param_name')
            if isinstance(node, ast.Call):
                if (  # pylint: disable=too-many-boolean-expressions
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "self"
                    and node.func.attr == "get"
                    and node.args
                    and isinstance(node.args[0], ast.Constant)
                    and isinstance(node.args[0].value, str)
                ):
                    params.add(node.args[0].value)

        return params


# ---------------------------------------------------------------------------
# Stage inference from validate_* roots and call graph
# ---------------------------------------------------------------------------

STAGE_ROOTS: Dict[str, List[str]] = {
    "common": ["validate_common"],
    "simulation": ["validate_simulation"],
    "pre_process": ["validate_pre_process"],
    "post_process": ["validate_post_process"],
}


def compute_method_stages(call_graph: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    """
    For each stage (simulation/pre_process/post_process/common), starting from
    validate_* roots, walk the call graph and record which methods belong to which stages.
    """
    method_stages: Dict[str, Set[str]] = defaultdict(set)

    def dfs(start: str, stage: str):
        stack = [start]
        visited: Set[str] = set()
        while stack:
            m = stack.pop()
            if m in visited:
                continue
            visited.add(m)
            method_stages[m].add(stage)
            for nxt in call_graph.get(m, ()):
                if nxt not in visited:
                    stack.append(nxt)

    for stage, roots in STAGE_ROOTS.items():
        for root in roots:
            dfs(root, stage)

    return method_stages


# ---------------------------------------------------------------------------
# Classification of messages for nicer grouping
# ---------------------------------------------------------------------------

def classify_message(msg: str) -> str:
    """
    Roughly classify rule messages for nicer grouping.

    Returns one of: "requirement", "incompatibility", "range", "other".
    """
    text = msg.lower()

    if (  # pylint: disable=too-many-boolean-expressions
        "not compatible" in text
        or "does not support" in text
        or "cannot be used" in text
        or "must not" in text
        or "is not supported" in text
        or "incompatible" in text
        or "untested" in text
    ):
        return "incompatibility"

    if (  # pylint: disable=too-many-boolean-expressions
        "requires" in text
        or "must be set if" in text
        or "must be specified" in text
        or "must be set with" in text
        or "can only be enabled if" in text
        or "must be set for" in text
    ):
        return "requirement"

    if (  # pylint: disable=too-many-boolean-expressions
        "must be between" in text
        or "must be positive" in text
        or "must be non-negative" in text
        or "must be greater than" in text
        or "must be less than" in text
        or "must be at least" in text
        or "must be <=" in text
        or "must be >=" in text
        or "must be odd" in text
        or "divisible by" in text
    ):
        return "range"

    return "other"


# Optional: nicer display names / categories (you can extend this)
FEATURE_META = {
    "igr": {"title": "Iterative Generalized Riemann (IGR)", "category": "solver"},
    "bubbles_euler": {"title": "Euler–Euler Bubble Model", "category": "bubbles"},
    "bubbles_lagrange": {"title": "Euler–Lagrange Bubble Model", "category": "bubbles"},
    "qbmm": {"title": "Quadrature-Based Moment Method (QBMM)", "category": "bubbles"},
    "polydisperse": {"title": "Polydisperse Bubble Dynamics", "category": "bubbles"},
    "mhd": {"title": "Magnetohydrodynamics (MHD)", "category": "physics"},
    "alt_soundspeed": {"title": "Alternative Sound Speed", "category": "physics"},
    "surface_tension": {"title": "Surface Tension Model", "category": "physics"},
    "hypoelasticity": {"title": "Hypoelasticity", "category": "physics"},
    "hyperelasticity": {"title": "Hyperelasticity", "category": "physics"},
    "relax": {"title": "Phase Change (Relaxation)", "category": "physics"},
    "viscous": {"title": "Viscosity", "category": "physics"},
    "acoustic_source": {"title": "Acoustic Sources", "category": "physics"},
    "ib": {"title": "Immersed Boundaries", "category": "geometry"},
    "cyl_coord": {"title": "Cylindrical Coordinates", "category": "geometry"},
    "weno_order": {"title": "WENO Order", "category": "numerics"},
    "muscl_order": {"title": "MUSCL Order", "category": "numerics"},
    "riemann_solver": {"title": "Riemann Solver", "category": "numerics"},
    "model_eqns": {"title": "Model Equations", "category": "fundamentals"},
    "num_fluids": {"title": "Number of Fluids", "category": "fundamentals"},
}


def feature_title(param: str) -> str:
    meta = FEATURE_META.get(param)
    if meta and "title" in meta:
        return meta["title"]
    return param


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
    """Get human-friendly model name"""
    models = {
        1: "π-γ (Compressible Euler)",
        2: "5-Equation (Multiphase)",
        3: "6-Equation (Phase Change)",
        4: "4-Equation (Single Component)"
    }
    return models.get(model_eqns, "Not specified")


def get_riemann_solver_name(solver: int | None) -> str:
    """Get Riemann solver name"""
    solvers = {
        1: "HLL",
        2: "HLLC",
        3: "Exact",
        4: "HLLD",
        5: "Lax-Friedrichs"
    }
    return solvers.get(solver, "Not specified")


def get_time_stepper_name(stepper: int | None) -> str:
    """Get time stepper name"""
    steppers = {
        1: "RK1 (Forward Euler)",
        2: "RK2",
        3: "RK3 (SSP)"
    }
    return steppers.get(stepper, "Not specified")


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
    lines.append(f"- [Model Equations (model_eqns = {summary['model_eqns']})](#-model-equations)")

    if summary['riemann_solver']:
        lines.append("- [Riemann Solvers](#️-riemann-solvers)")

    if summary['bubbles_euler'] or summary['bubbles_lagrange']:
        lines.append("- [Bubble Models](#-bubble-models)")

    if summary['mhd']:
        lines.append("- [MHD](#magnetohydrodynamics-mhd-mhd)")

    if summary['ib']:
        lines.append("- [Immersed Boundaries](#immersed-boundaries-ib)")

    if summary['viscous']:
        lines.append("- [Viscosity](#viscosity-viscous)")

    lines.append("\n</details>\n")
    return "\n".join(lines)


def generate_playbook() -> str:
    """Generate complete playbook from curated examples"""
    lines = []

    # Validate examples - will exit(1) if any are missing
    validate_playbook_examples()

    lines.append("## 🧩 Case Design Playbook\n")
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
    lines.append("## 🚀 Common Configuration Patterns\n")
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
    lines.append("## 📊 Feature Compatibility\n")
    lines.append("What works together:\n")

    for category, features in major_features.items():  # pylint: disable=too-many-nested-blocks
        lines.append(f"\n### {category}\n")

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

    # 3. Model Equations
    lines.append("## 🔢 Model Equations\n")
    lines.append("Choose your governing equations:\n")
    lines.append("")

    lines.append("<details>")
    lines.append("<summary><b>Model 1: π-γ (Compressible Euler)</b></summary>\n")
    lines.append("- **Use for:** Single-fluid compressible flow")
    lines.append("- **Value:** `model_eqns = 1`")
    lines.append("- **Note:** Cannot use `num_fluids`, bubbles, or certain WENO variants")
    lines.append("</details>\n")

    lines.append("<details>")
    lines.append("<summary><b>Model 2: 5-Equation (Most versatile)</b></summary>\n")
    lines.append("- **Use for:** Multiphase, bubbles, elastic materials, MHD")
    lines.append("- **Value:** `model_eqns = 2`")
    lines.append("- **Requirements:** Set `num_fluids`")
    lines.append("- **Compatible with:** Most physics models")
    lines.append("</details>\n")

    lines.append("<details>")
    lines.append("<summary><b>Model 3: 6-Equation (Phase change)</b></summary>\n")
    lines.append("- **Use for:** Phase change, cavitation")
    lines.append("- **Value:** `model_eqns = 3`")
    lines.append("- **Requirements:** `riemann_solver = 2` (HLLC), `avg_state = 2`, `wave_speeds = 1`")
    lines.append("- **Note:** Not compatible with bubbles or 3D cylindrical")
    lines.append("</details>\n")

    lines.append("<details>")
    lines.append("<summary><b>Model 4: 4-Equation (Single component)</b></summary>\n")
    lines.append("- **Use for:** Single-component flows with bubbles")
    lines.append("- **Value:** `model_eqns = 4`")
    lines.append("- **Requirements:** `num_fluids = 1`, set `rhoref` and `pref`")
    lines.append("</details>\n")

    # 4. Riemann Solvers (simplified)
    lines.append("## ⚙️ Riemann Solvers\n")
    lines.append("| Solver | `riemann_solver` | Best For | Requirements |")
    lines.append("|--------|-----------------|----------|-------------|")
    lines.append("| **HLL** | `1` | MHD, elastic materials | — |")
    lines.append("| **HLLC** | `2` | Bubbles, phase change, multiphase | `avg_state=2` for bubbles |")
    lines.append("| **Exact** | `3` | High accuracy (expensive) | — |")
    lines.append("| **HLLD** | `4` | MHD (advanced) | MHD only, no relativity |")
    lines.append("| **Lax-Friedrichs** | `5` | Robust fallback | Not with cylindrical+viscous |")
    lines.append("")

    # 5. Bubble Models (enhanced with collapsible)
    if "bubbles_euler" in by_param or "bubbles_lagrange" in by_param:
        lines.append("## 💧 Bubble Models\n")
        lines.append("")

        lines.append("<details>")
        lines.append("<summary><b>Euler-Euler (`bubbles_euler`)</b></summary>\n")
        lines.append("**Requirements:**")
        lines.append("- `model_eqns = 2` or `4`")
        lines.append("- `riemann_solver = 2` (HLLC)")
        lines.append("- `avg_state = 2`")
        lines.append("- Set `nb` (number of bins) ≥ 1\n")
        lines.append("**Extensions:**")
        lines.append("- `polydisperse = T`: Multiple bubble sizes (requires odd `nb > 1`)")
        lines.append("- `qbmm = T`: Quadrature method (requires `nnode = 4`)")
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

    # 6. Condensed Parameter Reference
    lines.append("## 📖 Quick Parameter Reference\n")
    lines.append("Key parameters and their constraints:\n")

    # Highlight only the most important parameters in collapsible sections
    important_params = {
        "MHD": "mhd",
        "Surface Tension": "surface_tension",
        "Viscosity": "viscous",
        "Number of Fluids": "num_fluids",
        "Cylindrical Coordinates": "cyl_coord",
        "Immersed Boundaries": "ib",
    }

    for title, param in important_params.items():
        if param not in by_param:
            continue

        rules_for_param = by_param[param]

        # Get key info
        requirements = []
        incompatibilities = []
        ranges = []

        for rule in rules_for_param:
            msg = rule.message
            # Skip IGR-related messages
            if "IGR" in msg:
                continue

            kind = classify_message(msg)
            if kind == "requirement":
                requirements.append(msg)
            elif kind == "incompatibility":
                incompatibilities.append(msg)
            elif kind == "range":
                ranges.append(msg)

        if not (requirements or incompatibilities or ranges):
            continue

        lines.append(f"\n<details>")
        lines.append(f"<summary><b>{title}</b> (`{param}`)</summary>\n")

        if requirements:
            lines.append("**Requirements:**")
            for req in requirements[:3]:
                lines.append(f"- {req}")
            lines.append("")

        if incompatibilities:
            lines.append("**Incompatibilities:**")
            for inc in incompatibilities[:3]:
                lines.append(f"- {inc}")
            lines.append("")

        if ranges:
            lines.append("**Valid values:**")
            for rng in ranges[:2]:
                lines.append(f"- {rng}")
            lines.append("")

        lines.append("</details>\n")

    # Add a footer with link to full validator
    lines.append("\n---\n")
    lines.append("💡 **Tip:** If you encounter a validation error, check the relevant section above or ")
    lines.append("review [`case_validator.py`](https://github.com/MFlowCode/MFC/blob/master/toolchain/mfc/case_validator.py) for complete validation logic.\n")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    src = CASE_VALIDATOR_PATH.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(CASE_VALIDATOR_PATH))

    analyzer = CaseValidatorAnalyzer()
    analyzer.visit(tree)

    # Infer stages per method from call graph
    method_stages = compute_method_stages(analyzer.call_graph)

    # Attach stages to rules
    for r in analyzer.rules:
        r.stages = method_stages.get(r.method, set())

    md = render_markdown(analyzer.rules)
    print(md)


if __name__ == "__main__":
    main()

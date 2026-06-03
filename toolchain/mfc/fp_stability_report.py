"""GitHub-output emitters for the FP-stability suite (step summary + annotations).

Pure formatting of the result dicts produced by the runners; the metric helpers
it uses (digit math) live in fp_stability_metrics.
"""

import math
import os

from .fp_stability_metrics import (
    MIN_SIG_BITS,
    VPREC_MANTISSA_BITS,
    _digits_left,
)


def _emit_github_annotations(results: list):
    """Emit GitHub annotations for FP cancellation sites.

    Only runs inside GitHub Actions (GITHUB_ACTIONS env var set). Annotations
    appear inline on the responsible source lines in the PR diff view.

    Up to 3 cancellation sites per case are emitted as ::notice:: so the diff
    highlights subtraction-cancellation hotspots from --check-cancellation. A site
    whose .fpp line sits inside a #:for/#:def expansion (tracked in
    cancellation_macro) is noted as possibly representing multiple instances.
    """
    if not os.environ.get("GITHUB_ACTIONS"):
        return
    for r in results:
        site_bits = r.get("cancellation_bits") or {}
        macro_sites = r.get("cancellation_macro") or {}
        for fname, lineno in r.get("cancellation_locs", [])[:3]:
            loc = f"file={fname},line={lineno}"
            title = f"FP cancellation [{r['name']}]"
            note = "catastrophic cancellation site"
            bits = site_bits.get((fname, lineno))
            if bits:
                note += f" — loses ≥ {bits / math.log2(10):.0f} of ~16 digits"
            macro = macro_sites.get((fname, lineno))
            if macro:
                note += f" — inside a {macro}-expanded line, may represent multiple instances"
            print(f"::notice {loc},title={title}::{note}", flush=True)
        n_cc = len(r.get("cancellation_locs", []))
        if n_cc > 3:
            print(f"::notice title=FP cancellation [{r['name']}]::{n_cc - 3} more cancellation site(s) not annotated inline; see the step summary", flush=True)


def _more_md(total: int, shown: int, noun: str) -> str:
    """Markdown bullet noting `total - shown` further items elided from a list,
    or '' when nothing was truncated."""
    if total <= shown:
        return ""
    return f"- …and {total - shown} more {noun}; see `fp-stability-logs/`"


def _emit_github_summary(results: list, n_samples: int):
    """Write a markdown results table to GITHUB_STEP_SUMMARY.

    Visible directly in the Actions run UI without downloading artifacts.
    Includes: pass/fail, max_dev, float proxy, VPREC sweep (failing levels),
    and catastrophic-cancellation source locations for any failing cases.
    """
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    n_pass = sum(1 for r in results if r["passed"])
    n_fail = len(results) - n_pass

    md = []
    md.append("## FP Stability Results\n")
    md.append(f"**{n_pass} passed, {n_fail} failed** — {n_samples} random-rounding samples per case\n")
    md.append(
        f"> **Coverage:** {len(results)} one-dimensional case(s) "
        f"({', '.join(r['name'] for r in results)}). A pass means stable in the code paths these "
        "cases exercise — not a guarantee for multi-D, viscous, MHD, IGR, or bubble-dynamics paths "
        "they do not reach.\n"
    )

    # Main results table — pass/fail is scale-free: bits retained vs a single floor
    md.append(f"_Pass = at least **{MIN_SIG_BITS} significant bits** retained under random rounding (scale-free; no per-case threshold)._\n")
    md.append("| Case | Status | bits retained | max\\_dev | Float proxy |")
    md.append("|------|:------:|:------:|--------:|--------:|")
    for r in results:
        status = "✅" if r["passed"] else "❌"
        bits = f"{r['sig_bits']:.1f}" if r.get("sig_bits") is not None else "—"
        fp = f"{r['float_proxy']:.2e}" if r["float_proxy"] is not None else "—"
        md.append(f"| `{r['name']}` | {status} | {bits} / {MIN_SIG_BITS} | {r['max_dev']:.2e} | {fp} |")
    md.append("")

    # Cancellation ORIGINS — where ill-conditioning actually arises, led with the
    # most severe (most bits lost).
    cases_with_cancel = [r for r in results if r.get("cancellation_locs")]
    if cases_with_cancel:
        md.append("### Catastrophic cancellation origins (ranked by digits lost)\n")
        md.append(
            "> Subtraction of nearly-equal values loses leading significant digits. A double carries "
            "~**16 significant digits** (53 bits); each entry shows how many that subtraction throws away "
            "(worst case, a lower bound). Losing ~8 digits halves your accuracy; losing ~13+ leaves only "
            "single-precision trust. Site *count* is not severity — one site losing many digits outweighs "
            "many mild ones.\n"
        )
        for r in cases_with_cancel:
            site_bits = r.get("cancellation_bits") or {}
            macro_sites = r.get("cancellation_macro") or {}
            sites = [{"where": f"{fname}:{lineno}", "bits": site_bits.get((fname, lineno), 0), "macro": macro_sites.get((fname, lineno))} for fname, lineno in r["cancellation_locs"]]
            ordered = sorted(sites, key=lambda e: (-e["bits"], e["where"]))
            if ordered:
                w = ordered[0]
                md.append(f"**`{r['name']}`** — {len(ordered)} site(s); worst loses ≥ {w['bits'] / math.log2(10):.0f} of ~16 digits\n")
            for e in ordered[:15]:
                lost = e["bits"] / math.log2(10)
                ambiguous = f" — _{e['macro']}-expanded, may represent multiple instances_" if e["macro"] else ""
                md.append(f"- **≥ {lost:.0f} digits lost** (~{_digits_left(e['bits']):.0f} of 16 left) — `{e['where']}`{ambiguous}")
            footer = _more_md(len(ordered), 15, "site(s)")
            if footer:
                md.append(footer)
            md.append("")

    # VPREC sweep — one column per mantissa-bit level showing the L∞ deviation at
    # that reduced precision (💥 crash = run diverged/failed, — = not measured).
    if any(r["vprec"] for r in results):
        _labels = {52: "52b", 23: "23b", 16: "16b", 10: "10b"}
        header = " | ".join(_labels[b] for b in VPREC_MANTISSA_BITS)
        sep = " | ".join(":---:" for _ in VPREC_MANTISSA_BITS)
        md.append("### VPREC precision sweep\n")
        md.append(f"| Case | {header} |")
        md.append(f"|------|{sep}|")
        for r in results:
            vmap = {b: d for b, d in r["vprec"]}
            cols = []
            for b in VPREC_MANTISSA_BITS:
                d = vmap.get(b)
                if d is None:
                    cols.append("—")
                elif d == float("inf"):
                    cols.append("💥 crash")
                else:
                    cols.append(f"{d:.2e}")
            md.append(f"| `{r['name']}` | {' | '.join(cols)} |")
        md.append("")

    # Float-max overflow sites
    cases_with_fmax = [r for r in results if r.get("float_max_locs")]
    if cases_with_fmax:
        md.append("### Float32 overflow sites (check\\_max\\_float)\n")
        for r in cases_with_fmax:
            md.append(f"**`{r['name']}`** — {len(r['float_max_locs'])} site(s)\n")
            for fname, lineno in r["float_max_locs"][:10]:
                md.append(f"- `{fname}:{lineno}`")
            footer = _more_md(len(r["float_max_locs"]), 10, "site(s)")
            if footer:
                md.append(footer)
            md.append("")

    with open(summary_path, "a") as f:
        f.write("\n".join(md) + "\n")

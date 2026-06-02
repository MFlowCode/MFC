"""GitHub-output emitters for the FP-stability suite (step summary + annotations).

Pure formatting of the result dicts produced by the runners; the metric helpers
it uses (statement resolution, source context, digit math) live in
fp_stability_metrics.
"""

import math
import os

from .fp_stability_metrics import (
    MIN_SIG_BITS,
    VPREC_MANTISSA_BITS,
    _digits_left,
    _get_source_context,
    _statement_at,
)


def _emit_github_annotations(results: list):
    """Emit GitHub annotations for FP hotspots.

    Only runs inside GitHub Actions (GITHUB_ACTIONS env var set). Annotations
    appear inline on the responsible source lines in the PR diff view.

    Up to 3 dd_line locations are emitted per case (minimal responsible lines
    from delta-debug).  Confirmed hotspots (suspect-only perturbation reproduced
    the instability) are ::warning::; unconfirmed ones are downgraded to
    ::notice:: so a suspect attribution is not presented as fact.  Up to 3
    cancellation sites per case are emitted as ::notice:: so the diff also
    highlights subtraction-cancellation hotspots from --check-cancellation.
    """
    if not os.environ.get("GITHUB_ACTIONS"):
        return
    for r in results:
        status = "FAIL" if not r["passed"] else "sensitivity"
        _sb = r.get("sig_bits")
        _sb_str = f"{_sb:.0f} bits retained (floor {MIN_SIG_BITS})" if _sb is not None else "n/a"
        dev_str = f"{_sb_str}, max_dev={r['max_dev']:.2e}"
        unconfirmed = r.get("dd_line_confirmed") is False

        for loc in r.get("dd_line_locs", [])[:3]:
            location = f"file={loc['path']},line={loc['start']}"
            if loc["end"] != loc["start"]:
                location += f",endLine={loc['end']}"
            note = dev_str
            if loc.get("share") is not None:
                note += f" — single-precision sensitivity: {loc['share'] * 100:.0f}% of float-proxy (where precision matters, not necessarily where cancellation originates)"
            if loc.get("cancellation"):
                note += " — also a catastrophic cancellation site"
            if loc.get("macro"):
                note += f" — {loc['macro']}-expanded line, may represent multiple instances"
            if unconfirmed:
                title = f"FP candidate (unconfirmed) [{r['name']}]"
                print(f"::notice {location},title={title}::{note}", flush=True)
            else:
                title = f"FP {status} [{r['name']}]"
                print(f"::warning {location},title={title}::{note}", flush=True)
        n_dd = len(r.get("dd_line_locs", []))
        if n_dd > 3:
            print(f"::notice title=FP hotspots [{r['name']}]::{n_dd - 3} more dd_line hotspot(s) not annotated inline; see the step summary", flush=True)

        for fname, lineno in r.get("cancellation_locs", [])[:3]:
            loc = f"file={fname},line={lineno}"
            title = f"FP cancellation [{r['name']}]"
            print(f"::notice {loc},title={title}::catastrophic cancellation site", flush=True)
        n_cc = len(r.get("cancellation_locs", []))
        if n_cc > 3:
            print(f"::notice title=FP cancellation [{r['name']}]::{n_cc - 3} more cancellation site(s) not annotated inline; see the step summary", flush=True)


def _more_md(total: int, shown: int, noun: str) -> str:
    """Markdown bullet noting `total - shown` further items elided from a list,
    or '' when nothing was truncated."""
    if total <= shown:
        return ""
    return f"- _…and {total - shown} more {noun}; see fp-stability-logs/_"


def _emit_github_summary(results: list, n_samples: int):
    """Write a markdown results table to GITHUB_STEP_SUMMARY.

    Visible directly in the Actions run UI without downloading artifacts.
    Includes: pass/fail, max_dev, float proxy, VPREC sweep (failing levels),
    and dd_line source locations for any failing cases.
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
    md.append("| Case | Status | bits retained | max\\_dev | Float proxy | MCA sig bits |")
    md.append("|------|:------:|:------:|--------:|--------:|:------:|")
    for r in results:
        status = "✅" if r["passed"] else "❌"
        bits = f"{r['sig_bits']:.1f}" if r.get("sig_bits") is not None else "—"
        fp = f"{r['float_proxy']:.2e}" if r["float_proxy"] is not None else "—"
        sb = str(r["mca_sigbits"]) if r.get("mca_sigbits") is not None else "—"
        md.append(f"| `{r['name']}` | {status} | {bits} / {MIN_SIG_BITS} | {r['max_dev']:.2e} | {fp} | {sb} |")
    md.append("")

    # Cancellation ORIGINS — where ill-conditioning actually arises, led with the
    # most severe (most bits lost). The numerically interesting signal; the
    # sensitivity list further down is dominated by the (benign) time integrator.
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
            # collapse continuation fragments to one entry per logical statement,
            # keeping the worst bits-lost seen on that statement
            stmts = {}  # (basename, stmt_start) -> {where, bits, text}
            for fname, lineno in r["cancellation_locs"]:
                stmt_start, _end, stmt_text = _statement_at(fname, lineno)
                key = (os.path.basename(fname), stmt_start)
                e = stmts.setdefault(key, {"where": f"{fname}:{stmt_start}", "bits": 0, "text": stmt_text})
                e["bits"] = max(e["bits"], site_bits.get((fname, lineno), 0))
            ordered = sorted(stmts.values(), key=lambda e: (-e["bits"], e["where"]))
            if ordered:
                w = ordered[0]
                md.append(f"**`{r['name']}`** — {len(stmts)} statement(s); worst loses ≥ {w['bits'] / math.log2(10):.0f} of ~16 digits\n")
            for e in ordered[:15]:
                lost = e["bits"] / math.log2(10)
                md.append(f"- **≥ {lost:.0f} digits lost** (~{_digits_left(e['bits']):.0f} of 16 left) — `{e['where']}`" + (f" — `{e['text']}`" if e["text"] else ""))
            footer = _more_md(len(ordered), 15, "statement(s)")
            if footer:
                md.append(footer)
            md.append("")

    # VPREC sweep — one column per bit level, ❌ where bits retained < floor
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

    # dd_line — single-precision SENSITIVITY (where precision most affects the
    # output). This is distinct from cancellation origin (reported separately):
    # the leader is typically the time integrator / final accumulation, because
    # perturbing the last write moves the output directly while upstream errors
    # get re-rounded there. Not a culprit-finder for ill-conditioning.
    cases_with_locs = [r for r in results if r["dd_line_locs"]]
    if cases_with_locs:
        md.append("<details>")
        md.append("<summary>Single-precision sensitivity (dd_line) — usually the time integrator; expand for details</summary>\n")
        md.append(
            "> Where reduced precision most moves the output — **typically the time integrator / "
            "final accumulation, which is expected and benign**. This is *not* where cancellation "
            "originates (that's the section above); it shows where precision matters most.\n"
        )
        _confirm_label = {True: "✅ confirmed", False: "⚠️ unconfirmed (suspect-only perturbation did not reproduce)", None: "— not checked"}
        for r in cases_with_locs:
            status = "❌ FAIL" if not r["passed"] else "✅ pass"
            md.append(f"**`{r['name']}`** ({status}) — attribution {_confirm_label[r.get('dd_line_confirmed')]}")
            md.append("_Ranked by the share of the single-precision deviation each line reproduces alone._\n")
            for loc in r["dd_line_locs"][:10]:
                rel_path, start, end = loc["path"], loc["start"], loc["end"]
                where = f"{rel_path}:{start}" if start == end else f"{rel_path}:{start}-{end}"
                tags = []
                if loc.get("share") is not None:
                    tags.append(f"**{loc['share'] * 100:.0f}%** of float-proxy")
                if loc.get("cancellation"):
                    tags.append("catastrophic cancellation")
                if loc.get("macro"):
                    tags.append(f"_{loc['macro']}-expanded, may represent multiple instances_")
                suffix = f" — {', '.join(tags)}" if tags else ""
                md.append(f"- `{where}`{suffix}")
                for inst in loc.get("instances", [])[:8]:
                    flag = " ⟵ flagrant" if inst is loc["instances"][0] and inst["dev"] > 0 else ""
                    md.append(f"  - instance #{inst['instance']} (`.f90:{inst['physline']}`, dev={inst['dev']:.2e}){flag}: `{inst['snippet']}`")
                snippet = _get_source_context(rel_path, start)
                if snippet:
                    md.append("  ```fortran")
                    for line in snippet.splitlines():
                        md.append(f"  {line}")
                    md.append("  ```")
            footer = _more_md(len(r["dd_line_locs"]), 10, "hotspot(s)")
            if footer:
                md.append(footer)
            md.append("")
        md.append("</details>\n")

    # dd_sym function names (collapsed, since less actionable than dd_line)
    cases_with_syms = [r for r in results if r["dd_sym_syms"]]
    if cases_with_syms:
        md.append("<details>")
        md.append("<summary>Responsible functions (dd_sym)</summary>\n")
        for r in cases_with_syms:
            md.append(f"\n**`{r['name']}`**\n")
            for sym in r["dd_sym_syms"]:
                md.append(f"- `{sym}`")
        md.append("\n</details>\n")

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

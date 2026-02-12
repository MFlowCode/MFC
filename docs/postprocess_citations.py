"""Post-process Doxygen HTML to add author names before citation numbers.

Transforms bare "[N]" citation links into "Author et al. [N]" style.
Parses docs/references.bib for author information and walks all HTML
files in the Doxygen output directory to insert author text.

Usage: python3 postprocess_citations.py <html_output_dir>
"""

import re
import sys
from pathlib import Path


def parse_bib_authors(bib_path: Path) -> dict[str, tuple[str, int, str | None]]:
    """Parse references.bib to extract first-author surname and author count.

    Returns:
        Dict mapping lowercase bib key ->
        (first_author_surname, author_count, second_author_surname | None)
    """
    text = bib_path.read_text(encoding="utf-8")
    entries: dict[str, tuple[str, int, str | None]] = {}

    # Find all bib entries: @type{Key,
    for m in re.finditer(r"@\w+\{(\w+),", text):
        key = m.group(1)
        # Find the author field for this entry
        start = m.end()
        # Look for author = {..} or author = "..."
        author_match = re.search(
            r"author\s*=\s*\{([^{}]*(?:\{[^}]*\}[^{}]*)*)\}",
            text[start:start + 2000],
        )
        if not author_match:
            continue

        author_str = author_match.group(1)
        # Split by " and " to get individual authors
        authors = re.split(r"\s+and\s+", author_str)
        count = len(authors)

        # Extract surname of first author
        first = authors[0].strip()
        surname = _extract_surname(first)

        # Also get second author surname for 2-author papers
        second_surname = None
        if count == 2:
            second_surname = _extract_surname(authors[1].strip())

        entries[key.lower()] = (surname, count, second_surname)

    return entries


def _extract_surname(author: str) -> str:
    """Extract surname from a BibTeX author name.

    Handles:
      - "G. Allaire" -> "Allaire"
      - "O. {Le M\\'etayer}" -> "Le Métayer"
      - "{Lord Rayleigh}" -> "Lord Rayleigh"
      - "B. {van Leer}" -> "van Leer"
      - "M. {Rodriguez Jr.}" -> "Rodriguez Jr."
    """
    # If the whole name is braced, return it (e.g., {Lord Rayleigh})
    if author.startswith("{") and author.endswith("}"):
        return _clean_latex(author[1:-1])

    # Check for braced surname at the end: "O. {Le Métayer}"
    brace_match = re.search(r"\{([^}]+)\}\s*$", author)
    if brace_match:
        return _clean_latex(brace_match.group(1))

    # Standard case: last word is surname
    parts = author.split()
    if parts:
        return _clean_latex(parts[-1])
    return author


def _clean_latex(s: str) -> str:
    """Remove LaTeX formatting from a string."""
    s = s.replace("\\'", "")  # accent commands
    s = s.replace('\\"', "")
    s = s.replace("\\", "")
    s = s.replace("{", "").replace("}", "")
    return s.strip()


# Pattern matching Doxygen citation links
CITE_RE = re.compile(
    r'<a class="el" href="citelist\.html#CITEREF_(\w+)">\[(\d+)\]</a>'
)


def process_html(html_path: Path, authors: dict) -> bool:
    """Process a single HTML file, inserting author names before citations.

    Returns True if the file was modified.
    """
    text = html_path.read_text(encoding="utf-8")
    if "citelist.html#CITEREF_" not in text:
        return False

    def replace_cite(match: re.Match) -> str:
        key = match.group(1)
        full_link = match.group(0)

        info = authors.get(key)
        if not info:
            return full_link

        surname, count, second_surname = info

        # Build author prefix
        if count == 1:
            prefix = surname
        elif count == 2 and second_surname:
            prefix = f"{surname} and {second_surname}"
        else:
            prefix = f"{surname} et al."

        # Dedup: only skip if this exact prefix was already inserted
        # (prevents double-insertion on re-runs of the script)
        preceding = text[max(0, match.start() - len(prefix) - 5):match.start()]
        if preceding.rstrip().endswith(prefix):
            return full_link

        return f"{prefix} {full_link}"

    new_text = CITE_RE.sub(replace_cite, text)
    if new_text != text:
        html_path.write_text(new_text, encoding="utf-8")
        return True
    return False


def check_bare_citations(html_dir: Path, authors: dict) -> list[str]:
    """Find citation links that lack an author-name prefix.

    Returns a list of warning strings for each bare citation found.
    """
    bare = []
    for html_file in sorted(html_dir.rglob("*.html")):
        if html_file.name == "citelist.html":
            continue
        text = html_file.read_text(encoding="utf-8")
        for m in CITE_RE.finditer(text):
            key = m.group(1)
            if key not in authors:
                bare.append(f"  {html_file.name}: [{m.group(2)}] (CITEREF_{key}) — no bib entry")
                continue
            # Check if an author name precedes the link
            before = text[max(0, m.start() - 60):m.start()].rstrip()
            surname, _count, second_surname = authors[key]
            has_prefix = (
                before.endswith(surname) or
                "et al." in before[-20:] or
                (second_surname and before.endswith(second_surname))
            )
            if not has_prefix:
                bare.append(f"  {html_file.name}: [{m.group(2)}] (CITEREF_{key}) — missing author prefix")
    return bare


def main():
    check_only = "--check" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--check"]

    if len(args) != 1:
        print(f"Usage: {sys.argv[0]} [--check] <html_output_dir>", file=sys.stderr)
        sys.exit(1)

    html_dir = Path(args[0])
    if not html_dir.is_dir():
        print(f"Error: {html_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    # Find references.bib relative to this script
    script_dir = Path(__file__).resolve().parent
    bib_path = script_dir / "references.bib"
    if not bib_path.exists():
        print(f"Error: {bib_path} not found", file=sys.stderr)
        sys.exit(1)

    authors = parse_bib_authors(bib_path)
    print(f"  Parsed {len(authors)} bib entries")

    if not check_only:
        modified = 0
        html_files = list(html_dir.rglob("*.html"))
        for html_file in html_files:
            if html_file.name == "citelist.html":
                continue
            if process_html(html_file, authors):
                modified += 1
        print(f"  Updated citations in {modified}/{len(html_files)} HTML files")

    # Verify all citations have author prefixes
    bare = check_bare_citations(html_dir, authors)
    if bare:
        print("Bare citations found (missing author prefix):")
        for b in bare:
            print(b)
        sys.exit(1)
    print("  All citations have author prefixes")


if __name__ == "__main__":
    main()

"""
IDE Configuration Module.

Automatically configures IDE settings (VS Code, etc.) for MFC development.
"""
# pylint: disable=import-outside-toplevel

import re
from pathlib import Path

from .common import MFC_ROOT_DIR

# Marker comments for the auto-generated section
_VSCODE_MARKER_BEGIN = "// MFC-SCHEMA-CONFIG-BEGIN (auto-generated, do not edit)"
_VSCODE_MARKER_END = "// MFC-SCHEMA-CONFIG-END"

# The MFC schema configuration to insert
# Matches common case file names - users get auto-completion for JSON/YAML case files
_VSCODE_MFC_CONFIG = '''\
    "json.schemas": [
        {
            "fileMatch": ["**/case.json", "**/input.json", "**/mfc-case.json", "**/mfc.json"],
            "url": "./toolchain/mfc-case-schema.json"
        }
    ],
    "yaml.schemas": {
        "./toolchain/mfc-case-schema.json": ["**/case.yaml", "**/case.yml", "**/input.yaml", "**/input.yml", "**/mfc-case.yaml", "**/mfc.yaml"]
    }'''


def ensure_vscode_settings() -> bool:
    """
    Ensure VS Code settings include MFC schema configuration.

    This is called on every mfc.sh invocation but is very lightweight:
    - Only reads/writes if the marker section is missing
    - Does not regenerate the schema (that's done via generate --json-schema)

    Returns:
        True if settings were updated, False if already configured
    """
    vscode_dir = Path(MFC_ROOT_DIR) / ".vscode"
    settings_path = vscode_dir / "settings.json"

    # Check if schema file exists (it should be committed to repo)
    schema_path = Path(MFC_ROOT_DIR) / "toolchain" / "mfc-case-schema.json"
    if not schema_path.exists():
        # Schema not generated yet - skip configuration
        return False

    # Build the marked config block
    marked_config = f"{_VSCODE_MARKER_BEGIN}\n{_VSCODE_MFC_CONFIG}\n    {_VSCODE_MARKER_END}"

    if settings_path.exists():
        content = settings_path.read_text()

        # Check if our markers already exist - if so, nothing to do
        if _VSCODE_MARKER_BEGIN in content:
            return False

        # Insert before the final closing brace
        last_brace = content.rfind('}')
        if last_brace != -1:
            # Check if we need a comma
            before_brace = content[:last_brace].rstrip()
            needs_comma = before_brace and not before_brace.endswith('{') and not before_brace.endswith(',')
            comma = ',' if needs_comma else ''
            new_content = (
                content[:last_brace].rstrip() +
                comma + '\n\n    ' +
                marked_config + '\n' +
                content[last_brace:]
            )
        else:
            # Malformed JSON, just append
            new_content = content + '\n' + marked_config
    else:
        # Ensure .vscode directory exists
        vscode_dir.mkdir(exist_ok=True)
        # Create new settings file with just our config
        new_content = f'{{\n    {marked_config}\n}}\n'

    settings_path.write_text(new_content)
    return True


def update_vscode_settings() -> None:
    """
    Force update VS Code settings with MFC schema configuration.

    Unlike ensure_vscode_settings(), this always updates the marked section,
    even if it already exists. Used by `generate --json-schema`.
    """
    from .printer import cons

    vscode_dir = Path(MFC_ROOT_DIR) / ".vscode"
    settings_path = vscode_dir / "settings.json"

    # Ensure .vscode directory exists
    vscode_dir.mkdir(exist_ok=True)

    # Build the marked config block
    marked_config = f"{_VSCODE_MARKER_BEGIN}\n{_VSCODE_MFC_CONFIG}\n    {_VSCODE_MARKER_END}"

    if settings_path.exists():
        content = settings_path.read_text()

        # Check if our markers already exist
        marker_pattern = re.compile(
            rf'{re.escape(_VSCODE_MARKER_BEGIN)}.*?{re.escape(_VSCODE_MARKER_END)}',
            re.DOTALL
        )

        if marker_pattern.search(content):
            # Replace existing marked section
            new_content = marker_pattern.sub(marked_config, content)
        else:
            # Insert before the final closing brace
            last_brace = content.rfind('}')
            if last_brace != -1:
                before_brace = content[:last_brace].rstrip()
                needs_comma = before_brace and not before_brace.endswith('{') and not before_brace.endswith(',')
                comma = ',' if needs_comma else ''
                new_content = (
                    content[:last_brace].rstrip() +
                    comma + '\n\n    ' +
                    marked_config + '\n' +
                    content[last_brace:]
                )
            else:
                new_content = content + '\n' + marked_config
    else:
        new_content = f'{{\n    {marked_config}\n}}\n'

    settings_path.write_text(new_content)
    cons.print(f"[green]Updated[/green] {settings_path}")

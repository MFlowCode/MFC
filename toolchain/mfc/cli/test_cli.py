"""
Smoke tests for cli/ modules.

Verifies that modules can be imported and basic functionality works.
"""
# pylint: disable=import-outside-toplevel

import unittest


class TestCliImports(unittest.TestCase):
    """Test that all CLI modules can be imported."""

    def test_schema_import(self):
        """Schema module should import and export expected classes."""
        from . import schema
        self.assertTrue(hasattr(schema, 'Command'))
        self.assertTrue(hasattr(schema, 'Argument'))
        self.assertTrue(hasattr(schema, 'Positional'))
        self.assertTrue(hasattr(schema, 'CLISchema'))

    def test_commands_import(self):
        """Commands module should import and have MFC_CLI_SCHEMA."""
        from . import commands
        self.assertTrue(hasattr(commands, 'MFC_CLI_SCHEMA'))
        self.assertIsNotNone(commands.MFC_CLI_SCHEMA)

    def test_argparse_gen_import(self):
        """Argparse generator should import."""
        from . import argparse_gen
        self.assertTrue(hasattr(argparse_gen, 'generate_parser'))

    def test_completion_gen_import(self):
        """Completion generator should import."""
        from . import completion_gen
        self.assertTrue(hasattr(completion_gen, 'generate_bash_completion'))
        self.assertTrue(hasattr(completion_gen, 'generate_zsh_completion'))

    def test_docs_gen_import(self):
        """Docs generator should import."""
        from . import docs_gen
        self.assertTrue(hasattr(docs_gen, 'generate_cli_reference'))


class TestCliSchema(unittest.TestCase):
    """Test CLI schema structure."""

    def test_cli_schema_has_commands(self):
        """MFC_CLI_SCHEMA should have commands defined."""
        from .commands import MFC_CLI_SCHEMA
        self.assertTrue(len(MFC_CLI_SCHEMA.commands) > 0)

    def test_cli_schema_has_description(self):
        """MFC_CLI_SCHEMA should have a description."""
        from .commands import MFC_CLI_SCHEMA
        self.assertIsNotNone(MFC_CLI_SCHEMA.description)
        self.assertIsInstance(MFC_CLI_SCHEMA.description, str)

    def test_commands_have_names(self):
        """Each command should have a name."""
        from .commands import MFC_CLI_SCHEMA
        for cmd in MFC_CLI_SCHEMA.commands:
            self.assertIsNotNone(cmd.name, f"Command missing name")
            self.assertTrue(len(cmd.name) > 0, f"Command has empty name")


class TestArgparseGenerator(unittest.TestCase):
    """Test argparse generator."""

    def test_generate_parser_returns_parser(self):
        """generate_parser should return a tuple with ArgumentParser."""
        import argparse
        from .argparse_gen import generate_parser
        from .commands import MFC_CLI_SCHEMA

        result = generate_parser(MFC_CLI_SCHEMA)
        # Returns (parser, subparsers_dict)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        parser, subparsers = result
        self.assertIsInstance(parser, argparse.ArgumentParser)
        self.assertIsInstance(subparsers, dict)

    def test_parser_has_subparsers(self):
        """Parser should have subparsers for each command."""
        from .argparse_gen import generate_parser
        from .commands import MFC_CLI_SCHEMA

        parser, subparsers = generate_parser(MFC_CLI_SCHEMA)
        # Should have subparsers for all commands
        self.assertTrue(len(subparsers) > 0)
        # Parser should not raise error when printing help
        try:
            parser.format_help()
        except Exception as e:
            self.fail(f"Parser help failed: {e}")


class TestCompletionGenerator(unittest.TestCase):
    """Test completion script generators."""

    def test_bash_completion_generates_output(self):
        """Bash completion should generate non-empty output."""
        from .completion_gen import generate_bash_completion
        from .commands import MFC_CLI_SCHEMA

        output = generate_bash_completion(MFC_CLI_SCHEMA)
        self.assertIsInstance(output, str)
        self.assertTrue(len(output) > 100)  # Should be substantial
        self.assertIn("complete", output.lower())  # Should contain bash complete

    def test_zsh_completion_generates_output(self):
        """Zsh completion should generate non-empty output."""
        from .completion_gen import generate_zsh_completion
        from .commands import MFC_CLI_SCHEMA

        output = generate_zsh_completion(MFC_CLI_SCHEMA)
        self.assertIsInstance(output, str)
        self.assertTrue(len(output) > 100)
        self.assertIn("compdef", output.lower())  # Should contain zsh compdef


class TestDocsGenerator(unittest.TestCase):
    """Test documentation generator."""

    def test_docs_generates_markdown(self):
        """Docs generator should produce markdown output."""
        from .docs_gen import generate_cli_reference
        from .commands import MFC_CLI_SCHEMA

        output = generate_cli_reference(MFC_CLI_SCHEMA)
        self.assertIsInstance(output, str)
        self.assertTrue(len(output) > 100)
        self.assertIn("#", output)  # Should contain markdown headers


if __name__ == "__main__":
    unittest.main()

"""
Tests for the viz module.

Covers step parsing, label formatting, format/timestep discovery,
data assembly (binary + silo, 1D/2D/3D), and 1D rendering.
Uses checked-in fixture data generated from minimal MFC runs.
"""
# pylint: disable=import-outside-toplevel

import os
import tempfile
import unittest


FIXTURES = os.path.join(os.path.dirname(__file__), 'fixtures')

# Fixture paths for each dimension + format
FIX_1D_BIN = os.path.join(FIXTURES, '1d_binary')
FIX_1D_SILO = os.path.join(FIXTURES, '1d_silo')
FIX_2D_BIN = os.path.join(FIXTURES, '2d_binary')
FIX_2D_SILO = os.path.join(FIXTURES, '2d_silo')
FIX_3D_BIN = os.path.join(FIXTURES, '3d_binary')
FIX_3D_SILO = os.path.join(FIXTURES, '3d_silo')


# ---------------------------------------------------------------------------
# Tests: _parse_steps
# ---------------------------------------------------------------------------

class TestParseSteps(unittest.TestCase):
    """Test _parse_steps() step-argument parsing."""

    def _parse(self, arg, available):
        from .viz import _parse_steps
        return _parse_steps(arg, available)

    def test_all_keyword(self):
        """'all' returns every available step."""
        self.assertEqual(self._parse('all', [0, 100, 200]), [0, 100, 200])

    def test_none_returns_all(self):
        """None returns every available step."""
        self.assertEqual(self._parse(None, [0, 100, 200]), [0, 100, 200])

    def test_last_keyword(self):
        """'last' returns only the final step."""
        self.assertEqual(self._parse('last', [0, 100, 200]), [200])

    def test_last_empty(self):
        """'last' with no available steps returns empty list."""
        self.assertEqual(self._parse('last', []), [])

    def test_single_int(self):
        """Single integer selects that step."""
        self.assertEqual(self._parse('100', [0, 100, 200]), [100])

    def test_single_int_missing(self):
        """Single integer not in available returns empty list."""
        self.assertEqual(self._parse('999', [0, 100, 200]), [])

    def test_range(self):
        """Range format 'start:end:stride' selects matching steps."""
        result = self._parse('0:500:200', [0, 100, 200, 300, 400, 500])
        self.assertEqual(result, [0, 200, 400])

    def test_range_no_stride(self):
        """Range format 'start:end' defaults to stride 1."""
        result = self._parse('0:2', [0, 1, 2, 3])
        self.assertEqual(result, [0, 1, 2])

    def test_invalid_value(self):
        """Non-numeric, non-keyword input raises MFCException."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException):
            self._parse('bogus', [0, 100])


# ---------------------------------------------------------------------------
# Tests: pretty_label
# ---------------------------------------------------------------------------

class TestPrettyLabel(unittest.TestCase):
    """Test pretty_label() LaTeX label generation."""

    def _label(self, varname):
        from .renderer import pretty_label
        return pretty_label(varname)

    def test_known_scalar(self):
        """Known scalars map to LaTeX."""
        self.assertEqual(self._label('pres'), r'$p$')
        self.assertEqual(self._label('rho'), r'$\rho$')

    def test_vel_indexed(self):
        """vel1/vel2/vel3 map to u/v/w."""
        self.assertEqual(self._label('vel1'), r'$u$')
        self.assertEqual(self._label('vel2'), r'$v$')
        self.assertEqual(self._label('vel3'), r'$w$')

    def test_alpha_indexed(self):
        """alpha<N> maps to LaTeX subscript."""
        self.assertIn('2', self._label('alpha2'))

    def test_unknown_passthrough(self):
        """Unknown variable names pass through unchanged."""
        self.assertEqual(self._label('my_custom_var'), 'my_custom_var')


# ---------------------------------------------------------------------------
# Tests: discover_format
# ---------------------------------------------------------------------------

class TestDiscoverFormat(unittest.TestCase):
    """Test discover_format() binary/silo detection."""

    def test_binary_detection(self):
        """Detects binary format from binary/ directory."""
        from .reader import discover_format
        self.assertEqual(discover_format(FIX_1D_BIN), 'binary')

    def test_silo_detection(self):
        """Detects silo format from silo_hdf5/ directory."""
        from .reader import discover_format
        self.assertEqual(discover_format(FIX_1D_SILO), 'silo')

    def test_missing_dir_raises(self):
        """Missing directories raise FileNotFoundError."""
        from .reader import discover_format
        d = tempfile.mkdtemp()
        try:
            with self.assertRaises(FileNotFoundError):
                discover_format(d)
        finally:
            os.rmdir(d)


# ---------------------------------------------------------------------------
# Tests: discover_timesteps
# ---------------------------------------------------------------------------

class TestDiscoverTimesteps(unittest.TestCase):
    """Test discover_timesteps() against fixture data."""

    def test_binary_1d(self):
        """Finds sorted timesteps from 1D binary fixture."""
        from .reader import discover_timesteps
        steps = discover_timesteps(FIX_1D_BIN, 'binary')
        self.assertEqual(steps, sorted(steps))
        self.assertIn(0, steps)
        self.assertGreater(len(steps), 1)

    def test_silo_1d(self):
        """Finds sorted timesteps from 1D silo fixture."""
        from .reader import discover_timesteps
        steps = discover_timesteps(FIX_1D_SILO, 'silo')
        self.assertEqual(steps, sorted(steps))
        self.assertIn(0, steps)
        self.assertGreater(len(steps), 1)


# ---------------------------------------------------------------------------
# Tests: binary read + assemble (1D, 2D, 3D)
# ---------------------------------------------------------------------------

class TestAssembleBinary1D(unittest.TestCase):
    """Test binary reader with 1D fixture data."""

    def test_ndim(self):
        """1D fixture assembles with ndim=1."""
        from .reader import assemble
        data = assemble(FIX_1D_BIN, 0, 'binary')
        self.assertEqual(data.ndim, 1)

    def test_grid_and_vars(self):
        """1D fixture has non-empty grid and expected variables."""
        from .reader import assemble
        data = assemble(FIX_1D_BIN, 0, 'binary')
        self.assertGreater(len(data.x_cc), 0)
        self.assertIn('pres', data.variables)
        self.assertIn('vel1', data.variables)
        self.assertEqual(data.variables['pres'].shape, data.x_cc.shape)

    def test_var_filter(self):
        """Passing var= loads only that variable."""
        from .reader import assemble
        data = assemble(FIX_1D_BIN, 0, 'binary', var='pres')
        self.assertIn('pres', data.variables)
        self.assertNotIn('vel1', data.variables)


class TestAssembleBinary2D(unittest.TestCase):
    """Test binary reader with 2D fixture data."""

    def test_ndim(self):
        """2D fixture assembles with ndim=2."""
        from .reader import assemble
        data = assemble(FIX_2D_BIN, 0, 'binary')
        self.assertEqual(data.ndim, 2)

    def test_grid_shape(self):
        """2D fixture has 2D variable arrays matching grid."""
        from .reader import assemble
        data = assemble(FIX_2D_BIN, 0, 'binary')
        self.assertGreater(len(data.x_cc), 0)
        self.assertGreater(len(data.y_cc), 0)
        pres = data.variables['pres']
        self.assertEqual(pres.shape, (len(data.x_cc), len(data.y_cc)))


class TestAssembleBinary3D(unittest.TestCase):
    """Test binary reader with 3D fixture data."""

    def test_ndim(self):
        """3D fixture assembles with ndim=3."""
        from .reader import assemble
        data = assemble(FIX_3D_BIN, 0, 'binary')
        self.assertEqual(data.ndim, 3)

    def test_grid_shape(self):
        """3D fixture has 3D variable arrays matching grid."""
        from .reader import assemble
        data = assemble(FIX_3D_BIN, 0, 'binary')
        self.assertGreater(len(data.x_cc), 0)
        self.assertGreater(len(data.y_cc), 0)
        self.assertGreater(len(data.z_cc), 0)
        pres = data.variables['pres']
        self.assertEqual(pres.shape,
                         (len(data.x_cc), len(data.y_cc), len(data.z_cc)))


# ---------------------------------------------------------------------------
# Tests: silo read + assemble (1D, 2D, 3D)
# ---------------------------------------------------------------------------

class TestAssembleSilo1D(unittest.TestCase):
    """Test silo reader with 1D fixture data."""

    def test_ndim(self):
        """1D silo fixture assembles with ndim=1."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_1D_SILO, 0)
        self.assertEqual(data.ndim, 1)

    def test_grid_and_vars(self):
        """1D silo fixture has non-empty grid and expected variables."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_1D_SILO, 0)
        self.assertGreater(len(data.x_cc), 0)
        self.assertIn('pres', data.variables)
        self.assertEqual(data.variables['pres'].shape, data.x_cc.shape)


class TestAssembleSilo2D(unittest.TestCase):
    """Test silo reader with 2D fixture data."""

    def test_ndim(self):
        """2D silo fixture assembles with ndim=2."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_2D_SILO, 0)
        self.assertEqual(data.ndim, 2)

    def test_grid_shape(self):
        """2D silo fixture has 2D variable arrays matching grid."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_2D_SILO, 0)
        pres = data.variables['pres']
        self.assertEqual(pres.shape, (len(data.x_cc), len(data.y_cc)))


class TestAssembleSilo3D(unittest.TestCase):
    """Test silo reader with 3D fixture data."""

    def test_ndim(self):
        """3D silo fixture assembles with ndim=3."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_3D_SILO, 0)
        self.assertEqual(data.ndim, 3)

    def test_grid_shape(self):
        """3D silo fixture has 3D variable arrays matching grid."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_3D_SILO, 0)
        pres = data.variables['pres']
        self.assertEqual(pres.shape,
                         (len(data.x_cc), len(data.y_cc), len(data.z_cc)))


# ---------------------------------------------------------------------------
# Tests: binary vs silo consistency
# ---------------------------------------------------------------------------

class TestBinarySiloConsistency(unittest.TestCase):
    """Verify binary and silo readers produce consistent results."""

    def test_1d_same_grid(self):
        """Binary and silo 1D fixtures have the same grid."""
        from .reader import assemble
        from .silo_reader import assemble_silo
        import numpy as np
        bin_data = assemble(FIX_1D_BIN, 0, 'binary')
        silo_data = assemble_silo(FIX_1D_SILO, 0)
        np.testing.assert_allclose(bin_data.x_cc, silo_data.x_cc, atol=1e-10)

    def test_1d_same_vars(self):
        """Binary and silo 1D fixtures have the same variable names."""
        from .reader import assemble
        from .silo_reader import assemble_silo
        bin_data = assemble(FIX_1D_BIN, 0, 'binary')
        silo_data = assemble_silo(FIX_1D_SILO, 0)
        self.assertEqual(sorted(bin_data.variables.keys()),
                         sorted(silo_data.variables.keys()))


# ---------------------------------------------------------------------------
# Tests: 1D rendering (requires matplotlib/imageio)
# ---------------------------------------------------------------------------

class TestRender1D(unittest.TestCase):
    """Smoke test: render 1D plots from fixture data."""

    def test_render_png(self):
        """Renders a single-variable PNG that is non-empty."""
        from .reader import assemble
        from .renderer import render_1d
        data = assemble(FIX_1D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_1d(data.x_cc, data.variables['pres'], 'pres', 0, out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)

    def test_render_tiled_png(self):
        """Tiled render of all variables produces a non-empty PNG."""
        from .reader import assemble
        from .renderer import render_1d_tiled
        data = assemble(FIX_1D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_1d_tiled(data.x_cc, data.variables, 0, out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)


class TestRender2D(unittest.TestCase):
    """Smoke test: render a 2D PNG from fixture data."""

    def test_render_2d_png(self):
        """Renders a 2D colormap PNG that is non-empty."""
        from .reader import assemble
        from .renderer import render_2d
        data = assemble(FIX_2D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_2d(data.x_cc, data.y_cc, data.variables['pres'],
                      'pres', 0, out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)


class TestRender3DSlice(unittest.TestCase):
    """Smoke test: render a 3D slice PNG from fixture data."""

    def test_render_3d_slice_png(self):
        """Renders a 3D midplane-slice PNG that is non-empty."""
        from .reader import assemble
        from .renderer import render_3d_slice
        data = assemble(FIX_3D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_3d_slice(data, 'pres', 0, out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)


if __name__ == "__main__":
    unittest.main()

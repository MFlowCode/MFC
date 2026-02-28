"""
Tests for the viz module.

Covers step parsing, label formatting, format/timestep discovery,
data assembly (binary + silo, 1D/2D/3D), and 1D rendering.
Uses checked-in fixture data generated from minimal MFC runs.
"""
# pylint: disable=import-outside-toplevel,protected-access

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

    def test_comma_list(self):
        """Comma-separated list selects the intersection with available steps."""
        result = self._parse('0,100,200,1000', [0, 100, 200, 300, 1000])
        self.assertEqual(result, [0, 100, 200, 1000])

    def test_comma_list_filters_unavailable(self):
        """Comma list silently drops steps not in available."""
        result = self._parse('0,999', [0, 100, 200])
        self.assertEqual(result, [0])

    def test_ellipsis_expansion(self):
        """Ellipsis infers stride and expands the range."""
        result = self._parse('0,100,200,...,1000',
                             list(range(0, 1001, 100)))
        self.assertEqual(result, list(range(0, 1001, 100)))

    def test_ellipsis_partial_available(self):
        """Ellipsis expansion filters to only available steps."""
        # only even-numbered hundreds available
        avail = [0, 200, 400, 600, 800, 1000]
        result = self._parse('0,100,...,1000', avail)
        self.assertEqual(result, [0, 200, 400, 600, 800, 1000])

    def test_ellipsis_requires_two_prefix_values(self):
        """Ellipsis with only one prefix value raises MFCException."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException):
            self._parse('0,...,1000', [0, 100, 1000])

    def test_ellipsis_must_be_second_to_last(self):
        """Ellipsis not in second-to-last position raises MFCException."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException):
            self._parse('0,100,...,500,1000', [0, 100, 500, 1000])

    def test_invalid_value(self):
        """Non-numeric, non-keyword input raises MFCException."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException):
            self._parse('bogus', [0, 100])

    def test_hyphen_range_raises_clean_error(self):
        """'0-100' (hyphen instead of colon) raises MFCException, not raw ValueError."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException):
            self._parse('0-100', [0, 100])


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

    def test_1d_same_values(self):
        """Binary and silo 1D fixtures have the same variable values."""
        import numpy as np
        from .reader import assemble
        from .silo_reader import assemble_silo
        bin_data = assemble(FIX_1D_BIN, 0, 'binary')
        silo_data = assemble_silo(FIX_1D_SILO, 0)
        common = sorted(set(bin_data.variables) & set(silo_data.variables))
        self.assertGreater(len(common), 0, "No common variables to compare")
        for vname in common:
            np.testing.assert_allclose(
                bin_data.variables[vname],
                silo_data.variables[vname],
                rtol=1e-5, atol=1e-10,
                err_msg=f"Variable '{vname}' differs between binary and silo",
            )


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


# ---------------------------------------------------------------------------
# Tests: _steps_hint
# ---------------------------------------------------------------------------

class TestStepsHint(unittest.TestCase):
    """Test _steps_hint() step preview for error messages."""

    def _hint(self, steps, n=8):
        from .viz import _steps_hint
        return _steps_hint(steps, n)

    def test_empty(self):
        """Empty steps returns 'none found'."""
        self.assertEqual(self._hint([]), "none found")

    def test_short_list_shows_all(self):
        """Short list shows all steps without truncation."""
        result = self._hint([0, 100, 200])
        self.assertIn('0', result)
        self.assertIn('200', result)
        self.assertNotIn('...', result)

    def test_long_list_truncated(self):
        """Long list includes count and truncation marker."""
        steps = list(range(0, 2000, 100))   # 20 steps
        result = self._hint(steps, n=8)
        self.assertIn('...', result)
        self.assertIn('[20 total]', result)
        self.assertIn('0', result)      # head present
        self.assertIn('1900', result)   # tail present


# ---------------------------------------------------------------------------
# Tests: _validate_cmap
# ---------------------------------------------------------------------------

class TestValidateCmap(unittest.TestCase):
    """Test _validate_cmap() colormap validation."""

    def _validate(self, name):
        from .viz import _validate_cmap
        _validate_cmap(name)

    def test_known_cmaps_pass(self):
        """Known colormaps do not raise."""
        for name in ('viridis', 'plasma', 'coolwarm', 'gray'):
            with self.subTest(name=name):
                self._validate(name)

    def test_unknown_cmap_raises(self):
        """Unknown colormap raises MFCException."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException):
            self._validate('notacolormap_xyz_1234')

    def test_typo_suggests_correct(self):
        """Typo in colormap name raises MFCException suggesting the correct spelling."""
        from mfc.common import MFCException
        with self.assertRaises(MFCException) as ctx:
            self._validate('virids')   # typo of viridis
        self.assertIn('viridis', str(ctx.exception))


# ---------------------------------------------------------------------------
# Tests: bounded TUI cache
# ---------------------------------------------------------------------------

class TestTuiCache(unittest.TestCase):
    """Test that the shared step cache respects CACHE_MAX."""

    def setUp(self):
        import mfc.viz._step_cache as cache_mod
        self._mod = cache_mod
        cache_mod.clear()

    def tearDown(self):
        self._mod.clear()

    def _read(self, step):
        return f"data_{step}"

    def test_cache_stores_entry(self):
        """Loaded step is stored in cache."""
        self._mod.load(0, self._read)
        self.assertIn(0, self._mod._cache)

    def test_cache_hit_avoids_reload(self):
        """Second load of same step does not call read_func again."""
        calls = [0]
        def counting(step):
            calls[0] += 1
            return step
        self._mod.load(5, counting)
        self._mod.load(5, counting)
        self.assertEqual(calls[0], 1)

    def test_cache_evicts_oldest_at_cap(self):
        """Oldest entry is evicted when CACHE_MAX is exceeded."""
        cap = self._mod.CACHE_MAX
        for i in range(cap + 3):
            self._mod.load(i, self._read)
        self.assertLessEqual(len(self._mod._cache), cap)
        self.assertNotIn(0, self._mod._cache)         # first evicted
        self.assertIn(cap + 2, self._mod._cache)      # most recent kept

    def test_seed_clears_and_populates(self):
        """seed() clears existing cache and pre-loads one entry."""
        self._mod.load(99, self._read)   # put something in first
        self._mod.seed(0, "preloaded")
        self.assertEqual(len(self._mod._cache), 1)
        self.assertEqual(self._mod._cache[0], "preloaded")


# ---------------------------------------------------------------------------
# Tests: log scale rendering (new feature smoke tests)
# ---------------------------------------------------------------------------

class TestRenderLogScale(unittest.TestCase):
    """Smoke test: log scale option produces valid PNG output."""

    def test_render_1d_log_scale(self):
        """render_1d with log_scale=True produces a non-empty PNG."""
        from .reader import assemble
        from .renderer import render_1d
        data = assemble(FIX_1D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_1d(data.x_cc, data.variables['pres'], 'pres', 0, out,
                      log_scale=True)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)

    def test_render_2d_log_scale(self):
        """render_2d with log_scale=True produces a non-empty PNG."""
        from .reader import assemble
        from .renderer import render_2d
        data = assemble(FIX_2D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_2d(data.x_cc, data.y_cc, data.variables['pres'],
                      'pres', 0, out, log_scale=True)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)


# ---------------------------------------------------------------------------
# Tests: multi-rank assembly (ghost-cell deduplication)
# ---------------------------------------------------------------------------

class TestMultiRankAssembly(unittest.TestCase):
    """Test assemble_from_proc_data with synthetic multi-processor data."""

    def _make_proc(self, x_cb, pres):
        """Build a minimal 1D ProcessorData from boundary coordinates."""
        import numpy as np
        from .reader import ProcessorData
        return ProcessorData(
            m=len(x_cb) - 1,
            n=0,
            p=0,
            x_cb=np.array(x_cb, dtype=np.float64),
            y_cb=np.array([0.0]),
            z_cb=np.array([0.0]),
            variables={'pres': np.array(pres, dtype=np.float64)},
        )

    def test_two_rank_1d_dedup(self):
        """Two processors with one overlapping ghost cell assemble correctly."""
        import numpy as np
        from .reader import assemble_from_proc_data
        # Domain: 4 cells with centers at 0.125, 0.375, 0.625, 0.875
        # Proc 0 sees cells 0-2 (center 0.625 is ghost from proc 1)
        # Proc 1 sees cells 1-3 (center 0.375 is ghost from proc 0)
        p0 = self._make_proc([0.00, 0.25, 0.50, 0.75],
                             [1.0, 2.0, 3.0])   # centers: 0.125, 0.375, 0.625
        p1 = self._make_proc([0.25, 0.50, 0.75, 1.00],
                             [2.0, 3.0, 4.0])   # centers: 0.375, 0.625, 0.875

        result = assemble_from_proc_data([(0, p0), (1, p1)])

        self.assertEqual(result.ndim, 1)
        self.assertEqual(len(result.x_cc), 4)
        np.testing.assert_allclose(result.x_cc, [0.125, 0.375, 0.625, 0.875])
        np.testing.assert_allclose(result.variables['pres'], [1.0, 2.0, 3.0, 4.0])

    def test_large_extent_dedup(self):
        """Deduplication works correctly for large-extent domains (>1e6)."""
        import numpy as np
        from .reader import assemble_from_proc_data
        # Scale up by 1e7: extent=1e7, decimals = ceil(-log10(1e7)) + 12 = 5
        scale = 1e7
        p0 = self._make_proc(
            [0.00 * scale, 0.25 * scale, 0.50 * scale, 0.75 * scale],
            [1.0, 2.0, 3.0],
        )
        p1 = self._make_proc(
            [0.25 * scale, 0.50 * scale, 0.75 * scale, 1.00 * scale],
            [2.0, 3.0, 4.0],
        )
        result = assemble_from_proc_data([(0, p0), (1, p1)])
        self.assertEqual(len(result.x_cc), 4)
        np.testing.assert_allclose(
            result.variables['pres'], [1.0, 2.0, 3.0, 4.0]
        )

    def test_very_large_extent_dedup_negative_decimals(self):
        """Deduplication works for extent ~1e13 where decimals becomes negative.

        At scale=1e13: extent = 1e13, decimals = ceil(-log10(1e13)) + 12 = -1.
        np.round with negative decimals rounds to the nearest 10^|d|, so
        np.round(x, -1) rounds to the nearest 10.  Cell widths of 2.5e12
        are >> 10, so distinct cell-centers must not be collapsed.
        """
        import numpy as np
        from .reader import assemble_from_proc_data
        scale = 1e13
        p0 = self._make_proc(
            [0.00 * scale, 0.25 * scale, 0.50 * scale, 0.75 * scale],
            [1.0, 2.0, 3.0],
        )
        p1 = self._make_proc(
            [0.25 * scale, 0.50 * scale, 0.75 * scale, 1.00 * scale],
            [2.0, 3.0, 4.0],
        )
        result = assemble_from_proc_data([(0, p0), (1, p1)])
        self.assertEqual(len(result.x_cc), 4)
        np.testing.assert_allclose(
            result.variables['pres'], [1.0, 2.0, 3.0, 4.0]
        )


# ---------------------------------------------------------------------------
# Tests: render_2d_tiled
# ---------------------------------------------------------------------------

class TestRender2DTiled(unittest.TestCase):
    """Smoke test: render_2d_tiled produces a valid PNG from 2D fixture data."""

    def test_render_2d_tiled_png(self):
        """Tiled render of all 2D variables produces a non-empty PNG."""
        from .reader import assemble
        from .renderer import render_2d_tiled
        data = assemble(FIX_2D_BIN, 0, 'binary')
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_2d_tiled(data, 0, out)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)


# ---------------------------------------------------------------------------
# Tests: render_3d_slice non-default axes and selectors
# ---------------------------------------------------------------------------

class TestRender3DSliceAxes(unittest.TestCase):
    """Test render_3d_slice with non-default slice axes and selectors."""

    def setUp(self):
        from .reader import assemble
        self._data = assemble(FIX_3D_BIN, 0, 'binary')

    def _render(self, **kwargs):
        from .renderer import render_3d_slice
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            out = f.name
        try:
            render_3d_slice(self._data, 'pres', 0, out, **kwargs)
            self.assertTrue(os.path.isfile(out))
            self.assertGreater(os.path.getsize(out), 0)
        finally:
            os.unlink(out)

    def test_x_axis_slice(self):
        """X-axis midplane slice produces a non-empty PNG."""
        self._render(slice_axis='x')

    def test_y_axis_slice(self):
        """Y-axis midplane slice produces a non-empty PNG."""
        self._render(slice_axis='y')

    def test_slice_by_index(self):
        """slice_index=0 selects first plane along default z axis."""
        self._render(slice_index=0)

    def test_slice_by_value(self):
        """slice_value selects the plane nearest the given coordinate."""
        z_mid = float(self._data.z_cc[len(self._data.z_cc) // 2])
        self._render(slice_value=z_mid)


# ---------------------------------------------------------------------------
# Tests: render_mp4
# ---------------------------------------------------------------------------

class TestRenderMp4(unittest.TestCase):
    """Smoke test: render_mp4 exercises frame rendering and returns a bool."""

    def _make_read_func(self, case_dir, fmt):
        from .reader import assemble
        def _read(step):
            return assemble(case_dir, step, fmt)
        return _read

    def test_mp4_1d_returns_bool(self):
        """render_mp4 with 1D data returns True or False without raising."""
        from .reader import discover_timesteps
        from .renderer import render_mp4
        steps = discover_timesteps(FIX_1D_BIN, 'binary')[:2]
        read_func = self._make_read_func(FIX_1D_BIN, 'binary')
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, 'test.mp4')
            result = render_mp4('pres', steps, out, fps=2, read_func=read_func)
            self.assertIsInstance(result, bool)

    def test_mp4_tiled_1d_returns_bool(self):
        """render_mp4 with tiled=True returns True or False without raising."""
        from .reader import discover_timesteps
        from .renderer import render_mp4
        steps = discover_timesteps(FIX_1D_BIN, 'binary')[:2]
        read_func = self._make_read_func(FIX_1D_BIN, 'binary')
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, 'test_tiled.mp4')
            result = render_mp4('pres', steps, out, fps=2,
                                read_func=read_func, tiled=True)
            self.assertIsInstance(result, bool)

    def test_mp4_no_read_func_raises(self):
        """render_mp4 with read_func=None raises ValueError."""
        from .renderer import render_mp4
        with self.assertRaises(ValueError):
            render_mp4('pres', [0], '/tmp/unused.mp4', read_func=None)

    def test_mp4_empty_steps_raises(self):
        """render_mp4 with empty steps raises ValueError."""
        from .renderer import render_mp4
        with self.assertRaises(ValueError):
            render_mp4('pres', [], '/tmp/unused.mp4',
                       read_func=lambda s: None)


# ---------------------------------------------------------------------------
# Tests: silo assemble_silo var_filter
# ---------------------------------------------------------------------------

class TestAssembleSiloVarFilter(unittest.TestCase):
    """Test assemble_silo with var= filter to cover the silo var_filter path."""

    def test_1d_var_filter_includes_only_requested(self):
        """Silo 1D: var='pres' loads pres and excludes vel1."""
        from .silo_reader import assemble_silo
        data = assemble_silo(FIX_1D_SILO, 0, var='pres')
        self.assertIn('pres', data.variables)
        self.assertNotIn('vel1', data.variables)

    def test_2d_var_filter_includes_only_requested(self):
        """Silo 2D: var='pres' loads pres and excludes other variables."""
        from .silo_reader import assemble_silo
        filtered = assemble_silo(FIX_2D_SILO, 0, var='pres')
        all_data = assemble_silo(FIX_2D_SILO, 0)
        self.assertIn('pres', filtered.variables)
        other_vars = [v for v in all_data.variables if v != 'pres']
        if other_vars:
            self.assertNotIn(other_vars[0], filtered.variables)


if __name__ == "__main__":
    unittest.main()

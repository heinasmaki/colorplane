"""Numerical and command-line checks for the optional plotting tools."""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lattice_visualization import (crop_curve, interaction_kernel, load_lattice,
                                   load_probability_map, neighbor_probabilities, settings)


def pair_reference(lattice, spacing, kernel, normalization="full"):
    # Enumerate pairs directly; do not use the production stencil or FFT.
    radius = math.ceil(2 + 1 / spacing)

    def weight(dy, dx):
        r = math.hypot(dx, dy) - 1 / spacing
        value = (max(0, 1-abs(r)) if kernel == "hat" else
                 .25*(1+math.cos(math.pi*r/2)) if abs(r) < 2 else 0) / spacing
        return value if value > 1e-14 else 0

    norm = math.fsum(weight(dy, dx) for dy in range(-radius, radius+1)
                     for dx in range(-radius, radius+1))
    result = np.zeros(lattice.shape)
    for y, x in np.ndindex(lattice.shape):
        numerator = denominator = 0.
        for v, u in np.ndindex(lattice.shape):
            w = weight(v-y, u-x)
            denominator += w
            if lattice[y, x] == lattice[v, u]:
                numerator += w
        denominator = norm if normalization == "full" else denominator
        result[y, x] = numerator / denominator if denominator else np.nan
    return result


class VisualizationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="colorplane-plot-")
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_fft_matches_pair_enumeration_at_free_boundaries(self):
        lattice = np.random.default_rng(42).integers(1, 4, size=(9, 13))
        for kernel in ("cosine", "hat"):
            for normalization in ("full", "available"):
                actual = neighbor_probabilities(lattice, .2, kernel, normalization)
                expected = pair_reference(lattice, .2, kernel, normalization)
                np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-14)

    def test_uniform_interior_and_boundary_normalizations(self):
        lattice = np.ones((25, 27), dtype=int)
        full = neighbor_probabilities(lattice, .2)
        available = neighbor_probabilities(lattice, .2, normalization="available")
        self.assertLess(full[0, 0], .5)
        self.assertAlmostEqual(full[12, 13], 1, places=14)
        np.testing.assert_allclose(available, 1, atol=1e-14)
        tiny = neighbor_probabilities(np.ones((1, 1)), .2, normalization="available")
        self.assertTrue(np.isnan(tiny[0, 0]))

    def test_fft_does_not_wrap_and_supports_single_row(self):
        lattice = np.array([[1, 2, 3, 4, 5, 6, 7, 1]])
        np.testing.assert_allclose(neighbor_probabilities(lattice, .3),
                                   pair_reference(lattice, .3, "cosine"), atol=1e-14)

    def test_crop_curve_matches_direct_slices(self):
        values = np.arange(35, dtype=float).reshape(5, 7) / 35
        values[1, 1] = np.nan
        layers, means = crop_curve(values)
        for layer, mean in zip(layers, means):
            cropped = values[layer:5-layer, layer:7-layer]
            self.assertAlmostEqual(mean, np.nanmean(cropped), places=14)
        layers, means = crop_curve(np.full((2, 3), np.nan))
        self.assertEqual(layers.tolist(), [0])
        self.assertTrue(np.isnan(means[0]))

    def test_readers_and_metadata_overrides(self):
        path = self.directory / "run.lat"
        path.write_text("1 2 255\n")
        self.assertEqual(load_lattice(path).shape, (1, 3))
        path.with_suffix(".meta").write_text("spacing=2D-1\nkernel=hat\nq=255\n")
        self.assertEqual(settings(path), (.2, "hat", 255))
        self.assertEqual(settings(path, .1, "cosine", 3), (.1, "cosine", 3))
        path.with_suffix(".meta").unlink()
        path.with_suffix(".sum").write_text("Lattice spacing: 0.01250\n")
        self.assertEqual(settings(path)[0], .0125)
        map_path = path.with_suffix(".ngb")
        map_path.write_text("0D0 5D-1 1D0\n")
        np.testing.assert_array_equal(load_probability_map(map_path, (1, 3)), [[0, .5, 1]])
        with self.assertRaises(ValueError):
            load_probability_map(map_path, (3, 1))
        for text in ("0 1\n", "1.5 2\n", "1 NaN\n", "1 256\n", "1 2\n3\n"):
            path.write_text(text)
            with self.assertRaises(ValueError):
                load_lattice(path)

    def test_fortran_local_energy_map(self):
        probe = ROOT / "build" / "tests" / ("energy_probe.exe" if os.name == "nt" else "energy_probe")
        if not probe.exists():
            self.skipTest("Run tests/test_reliability.py first to build the optional Fortran probe")
        lattice = np.random.default_rng(37).integers(1, 4, size=(7, 9))
        path = self.directory / "run.lat"
        np.savetxt(path, lattice, fmt="%d")
        for kernel in ("cosine", "hat"):
            proc = subprocess.run([str(probe)], input=f"3 9 7 .2\n{kernel}\n{path}\n",
                                  capture_output=True, text=True, check=True)
            # Probe prints a summary first and then all local energies.
            actual = np.array([[float(v) for v in line.split()] for line in proc.stdout.splitlines()[1:]])
            np.testing.assert_allclose(neighbor_probabilities(lattice, .2, kernel), actual,
                                       rtol=1e-13, atol=1e-14)

    def run_script(self, script, *args, ok=True):
        proc = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
                              capture_output=True, text=True, timeout=120)
        self.assertEqual(proc.returncode == 0, ok, proc.stdout + proc.stderr)
        return proc

    def test_headless_clis_and_explicit_map_precedence(self):
        path = self.directory / "run.lat"
        np.savetxt(path, np.array([[1, 2, 1], [2, 1, 2]]), fmt="%d")
        png = self.directory / "lattice.png"
        self.run_script("plot_lattice.py", path, "--output", png)
        self.assertEqual(png.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
        map_path = path.with_suffix(".ngb")
        np.savetxt(map_path, np.full((2, 3), .25))
        exported = self.directory / "displayed.txt"
        self.run_script("plot_neighbors.py", path, "--output", png, "--data-output", exported)
        np.testing.assert_array_equal(np.loadtxt(exported), .25*np.ones((2, 3)))
        path.with_suffix(".meta").write_text("spacing=.2\nkernel=hat\nq=2\n")
        self.run_script("plot_neighbors.py", path, "--recompute", "--output", png,
                        "--data-output", exported)
        np.testing.assert_allclose(np.loadtxt(exported), neighbor_probabilities(load_lattice(path), .2, "hat"))
        self.run_script("plot_neighbors.py", path, "--normalization", "available",
                        "--output", png, "--data-output", exported)
        self.assertTrue(np.isnan(np.loadtxt(exported)).all())
        self.run_script("plot_neighbors.py", path, "--data-output", path, ok=False)
        self.run_script("plot_neighbors.py", path, "--map", map_path, "--normalization", "available", ok=False)
        self.run_script("plot_lattice.py", path, "--colors", 1, ok=False)


if __name__ == "__main__":
    unittest.main(verbosity=2)

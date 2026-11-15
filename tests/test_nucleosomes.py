"""Nucleosome interpolation (the class itself, independent of the pipeline)."""
import numpy as np

from multimm.nucleosome_interpolation import NucleosomeInterpolation


def line(n=6):
    return np.column_stack([np.arange(n) * 0.2, np.zeros(n), np.zeros(n)])


def test_max_nucleosomes_from_radius():
    assert NucleosomeInterpolation(line(), np.ones(6), nuc_radius=0.1).max_nucs_per_bead == 5
    assert NucleosomeInterpolation(line(), np.ones(6), nuc_radius=0.25).max_nucs_per_bead == 2


def test_no_signal_means_no_nucleosomes():
    V = line()
    out = NucleosomeInterpolation(V, np.zeros(6)).interpolate_structure_with_nucleosomes()
    np.testing.assert_allclose(out, V)  # beads only, nothing added


def test_signal_adds_helix_points():
    V = line()
    bw = np.array([0.0, 1.0, 5.0, 20.0, 50.0, 100.0])
    out = NucleosomeInterpolation(V, bw, points_per_nuc=20).interpolate_structure_with_nucleosomes()
    assert len(out) > len(V)
    assert (len(out) - len(V)) % 20 == 0  # whole nucleosomes of 20 points each

"""Pure helpers from utils.py."""
import numpy as np
import pytest

from multimm import utils as u


def test_min_max_normalize():
    out = u.min_max_normalize(np.array([[1.0, 3.0], [5.0, np.nan]]))
    assert out.min() == 0 and out.max() == 1


def test_compute_averages():
    assert np.allclose(u.compute_averages(np.arange(8.0), 4), [0.5, 2.5, 4.5, 6.5])


def test_rescale_and_downsample():
    m = np.arange(100.0).reshape(10, 10)
    assert u.rescale_matrix(m, 5).shape == (5, 5)
    assert u.mean_downsample(np.random.rand(40, 3), 10).shape == (10, 3)


def test_remove_zero_rows_and_columns():
    m = np.ones((4, 4))
    m[1] = 0
    m[:, 1] = 0
    out = u.remove_zero_rows_and_columns(m)
    assert out[0].shape == (3, 3)


def test_standarize():
    s = u.standarize(np.array([1.0, 2.0, 3.0, 4.0]))
    assert s.mean() == pytest.approx(0) and s.std() == pytest.approx(1)


def test_chrom_sort_key():
    assert sorted(["chrX", "chr10", "chr2", "chrY"], key=u.chrom_sort_key) == ["chr2", "chr10", "chrX", "chrY"]


def test_discretize_array():
    out = u.discretize_array(np.array([0.1, 0.3, 0.5, 0.7, 0.9]), [0.2, 0.4, 0.6, 0.8])
    assert len(out) == 5 and out[0] != out[-1]


def test_shuffle_blocks_preserves_content():
    a = np.array([1, 1, 2, 2, 2, 3])
    assert sorted(u.shuffle_blocks(a)) == sorted(a)


def test_generate_arrays():
    ms, ns = u.generate_arrays(20, 100)[:2]
    assert len(ms) == len(ns) == 20 and (ns >= ms).all()


def test_compartment_helpers():
    pc1 = np.array([-1.0, -0.5, 0.5, 1.0])
    dens = np.array([0.1, 0.2, 0.9, 1.0])  # dense = B
    aligned = u.align_pc1_sign(pc1, dens)
    assert aligned[dens.argmax()] < 0 < aligned[dens.argmin()]
    assert set(u.discretize_compartments(aligned)) <= {-1, 0, 1}


def test_pc1_and_oe():
    n = 40
    i = np.arange(n)
    m = 1 / (1 + np.abs(i[:, None] - i[None]))
    assert u.oe_matrix(m).shape == m.shape
    assert len(u.hic_pc1(m, already_oe=False)[0] if isinstance(u.hic_pc1(m), tuple) else u.hic_pc1(m)) == n


def test_integers_to_hex_colors():
    c = u.integers_to_hex_colors(1, 4)
    assert len(c) == 4 and all(x.startswith("#") for x in c)


def test_model_distance_heatmap():
    X = np.random.rand(12, 3)
    assert u.model_distance_heatmap(X).shape == (12, 12)

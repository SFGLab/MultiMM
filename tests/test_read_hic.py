"""Hi-C matrix preprocessing helpers."""
import numpy as np
import pytest

from multimm import read_hic as rh


@pytest.fixture()
def mat():
    rng = np.random.default_rng(0)
    a = rng.random((40, 40))
    return (a + a.T) / 2


def test_symmetrize_and_clean(mat):
    m = mat.copy()
    m[3, 5] = np.nan
    m[0, 1] = -1
    out = rh.symmetrize_and_clean(m)
    assert np.allclose(out, out.T)
    assert np.isfinite(out).all() and (out >= 0).all()
    assert np.allclose(np.diag(out), 0)


@pytest.mark.parametrize("n_out", [10, 40, 80])
def test_pool_to_n_beads_shape(mat, n_out):
    assert rh.pool_to_n_beads(mat, n_out).shape == (n_out, n_out)


def test_pool_weights_rows_sum_to_one():
    w = rh._pool_weights(40, 10)
    assert w.shape[0] == 10 or w.shape[1] == 10


def test_handle_missing_bins_fills_short_gap(mat):
    m = mat.copy()
    m[10, :] = 0
    m[:, 10] = 0
    out = rh.handle_missing_bins(m, max_gap=3)
    assert out.shape == m.shape
    assert out[10].sum() > 0


def test_oe_enrichment_is_finite(mat):
    out = rh.oe_enrichment_matrix(mat)
    assert out.shape == mat.shape and np.isfinite(out).all()


def test_diagonal_normalize_keeps_shape(mat):
    assert rh.diagonal_normalize(mat).shape == mat.shape


def test_denoise_bounds(mat):
    out = rh.denoise_contact_matrix(mat)
    assert out.shape == mat.shape and np.isfinite(out).all()


def test_preprocess_shape_and_range(mat):
    out = rh.preprocess_hic_matrix(mat, 20)
    assert out.shape == (20, 20)
    assert np.isfinite(out).all() and out.min() >= 0


def test_detect_format(tmp_path):
    for name, fmt in [("a.hic", "hic"), ("a.cool", "cool"), ("a.mcool", "mcool")]:
        p = tmp_path / name
        p.write_text("x")
        assert rh._detect_format(str(p)) == fmt
    bad = tmp_path / "a.txt"
    bad.write_text("x")
    with pytest.raises(ValueError):
        rh._detect_format(str(bad))
    with pytest.raises(FileNotFoundError):
        rh._detect_format(str(tmp_path / "missing.hic"))


def test_normalise_chrom():
    assert rh._normalise_chrom("1", ["chr1", "chr2"]) == "chr1"
    assert rh._normalise_chrom("chr2", ["1", "2"]) == "2"

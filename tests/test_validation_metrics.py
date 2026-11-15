"""Hi-C validation metrics on synthetic matrices."""
import numpy as np
import pytest

from multimm import validation as v


@pytest.fixture()
def plaid():
    # two-block plaid map with distance decay
    n = 60
    i = np.arange(n)
    decay = 1.0 / (1 + np.abs(i[:, None] - i[None, :]))
    blocks = (np.arange(n) // 15 % 2)
    plaid = 1 + 0.5 * (blocks[:, None] == blocks[None, :])
    return decay * plaid


def test_diagonal_decay_decreases(plaid):
    p = v.diagonal_decay_profile(plaid)
    assert len(p) == plaid.shape[0] - 1 and p[0] > p[-1]
    assert len(v.diagonal_decay_profile(plaid, max_diag=5)) == 5


def test_insulation_shape(plaid):
    assert v.insulation_score(plaid, window=5).shape[0] == plaid.shape[0]


def test_pc1_of_oe_shape(plaid):
    assert v.pc1_of_oe(plaid).shape == (plaid.shape[0],)


def test_inverse_contact_matrix_monotone():
    d = np.array([[0, 1.0], [1.0, 0]])
    m = v.inverse_contact_matrix(d)
    assert m.shape == d.shape and np.isfinite(m).all()


def test_pearson_spearman_perfect():
    a = np.arange(20, dtype=float)
    assert v._pearson(a, 2 * a + 1)[0] == pytest.approx(1.0)
    assert v._spearman(a, a ** 3)[0] == pytest.approx(1.0)


def test_normalisers_range():
    y = np.array([3.0, 5.0, 9.0])
    assert v._normalize01(y).min() == 0 and v._normalize01(y).max() == 1
    assert abs(v._normalize_signed(np.array([-2.0, 1.0, 4.0]))).max() <= 1 + 1e-9


def test_upper_tri_count():
    assert v._upper_tri(np.ones((5, 5))).size == 10


def test_similarity_metrics_identical_vs_noise(plaid):
    rng = np.random.default_rng(1)
    noise = rng.random(plaid.shape)
    assert v._ssim(plaid, plaid) == pytest.approx(1.0, abs=1e-6)
    assert v._ssim(plaid, plaid) > v._ssim(plaid, noise)
    assert v._gmsd(plaid, plaid) == pytest.approx(0.0, abs=1e-6)
    assert v._gmsd(plaid, plaid) < v._gmsd(plaid, noise)
    assert v._mutual_information(plaid, plaid) > v._mutual_information(plaid, noise)


def test_pool_and_smooth(plaid):
    assert v._pool_matrix(plaid, 20).shape == (20, 20)
    assert v._smooth_matrix(plaid).shape == plaid.shape
    assert v._smooth_signal(np.arange(10.0), 1.0).shape == (10,)


def test_coords_to_contact_matrices():
    xyz = np.random.default_rng(2).random((30, 3)).astype(np.float32)
    inv = v._coords_to_inv_contact_f32(xyz)
    assert inv.shape == (30, 30)

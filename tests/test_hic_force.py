"""Hi-C force: kernels, restraint construction and calibration."""
import numpy as np
import pytest

from multimm.hic_force import (
    VALID_BOLTZMANN_KERNELS,
    _boltzmann_r_target,
    auto_contact_scale,
    build_boltzmann_force,
    get_boltzmann_p_func,
)

RC, R_MIN, ALPHA = 1.0, 0.1, 4.0


@pytest.mark.parametrize("kernel", VALID_BOLTZMANN_KERNELS)
def test_target_distance_inverts_probability(kernel):
    c = np.linspace(0.1, 0.9, 9)
    r = _boltzmann_r_target(c, RC, ALPHA, R_MIN, kernel=kernel)
    inside = (r > R_MIN + 1e-9) & (r < RC - 1e-9)  # clipped values cannot invert
    assert inside.any()
    p = get_boltzmann_p_func(RC, ALPHA, R_MIN, kernel)(r[inside])
    np.testing.assert_allclose(p, c[inside], rtol=1e-6)


@pytest.mark.parametrize("kernel", VALID_BOLTZMANN_KERNELS)
def test_target_distance_is_clipped_and_monotone(kernel):
    c = np.array([1e-9, 0.01, 0.1, 0.5, 1.0])
    r = _boltzmann_r_target(c, RC, ALPHA, R_MIN, kernel=kernel)
    assert r.min() >= R_MIN - 1e-12 and r.max() <= RC + 1e-12
    assert np.all(np.diff(r) <= 1e-12)  # stronger contact -> shorter distance


def test_unknown_kernel_rejected():
    with pytest.raises(ValueError, match="KERNEL"):
        _boltzmann_r_target(np.array([0.5]), RC, ALPHA, R_MIN, kernel="nope")


def test_exponential_formula():
    r = _boltzmann_r_target(np.array([0.5]), RC, ALPHA, R_MIN, kernel="exponential")
    assert r[0] == pytest.approx(R_MIN - (RC / ALPHA) * np.log(0.5))


def test_capped_pairs_get_one_sided_floor():
    # one strong pair (real distance) and one very weak pair (capped at rc)
    C = np.zeros((4, 4))
    C[0, 1] = C[1, 0] = 1.0
    C[0, 3] = C[3, 0] = 1e-12
    force = build_boltzmann_force(C, rc=RC, alpha=ALPHA, k_scale=40.0, r_min=R_MIN)
    assert force.getNumBonds() == 2
    pars = {tuple(force.getBondParameters(i)[:2]): force.getBondParameters(i)[2] for i in range(2)}
    c, r_lo, r_hi = pars[(0, 1)]
    assert r_lo == pytest.approx(r_hi) and r_lo == pytest.approx(R_MIN)  # two-sided well
    c, r_lo, r_hi = pars[(0, 3)]
    assert r_lo == pytest.approx(RC) and r_hi == pytest.approx(10 * RC)  # floor only


def test_matrix_normalised_to_unit_maximum():
    C = np.zeros((3, 3))
    C[0, 2] = C[2, 0] = 50.0
    force = build_boltzmann_force(C, rc=RC, alpha=ALPHA, r_min=R_MIN)
    assert force.getBondParameters(0)[2][0] == pytest.approx(1.0)


def test_auto_contact_scale_is_requested_percentile():
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(500, 3))
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)[np.triu_indices(500, 1)]
    assert auto_contact_scale(xyz, percentile=50.0) == pytest.approx(np.median(d), rel=0.05)

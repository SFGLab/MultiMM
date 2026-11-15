"""B-lamina profiles: maximum mid-shell, minima at the walls."""
from types import SimpleNamespace

import numpy as np
import openmm as mm
import pytest

from multimm.model import MultiMM

R1, R2, B = 2.0, 6.0, 400.0


def _energies(mode, radii):
    # stub self: only what add_Blamina_interaction touches
    system = mm.System()
    for _ in radii:
        system.addParticle(1.0)
    stub = SimpleNamespace(
        args=SimpleNamespace(BLAMINA_FORCE_TYPE=mode, IBL_SCALE=B),
        radius1=R1, radius2=R2, mass_center=np.zeros(3),
        system=system, Cs=[-1] * len(radii),
        _register_force=lambda f, n: 0,
    )
    MultiMM.add_Blamina_interaction(stub)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions([mm.Vec3(r, 0, 0) for r in radii])
    # one particle at a time would need a context each; use a probe sweep
    return ctx.getState(getEnergy=True).getPotentialEnergy()._value


def _profile(mode, radii):
    return np.array([_energies(mode, [r]) for r in radii])


@pytest.mark.parametrize("mode", ["sin", "gaussian_shell", "harmonic_shell", "logistic_shell"])
def test_profile_shape(mode):
    r = np.linspace(0.5, 8.0, 76)
    e = _profile(mode, r)
    mid = 0.5 * (R1 + R2)
    # maximum inside the shell, near the middle
    inside = (r > R1) & (r < R2)
    assert abs(r[inside][np.argmax(e[inside])] - mid) < 0.5
    # minima near the walls
    left, right = r[(r > R1) & (r < mid)], r[(r > mid) & (r < R2)]
    assert abs(left[np.argmin(e[(r > R1) & (r < mid)])] - R1) < 0.8
    assert abs(right[np.argmin(e[(r > mid) & (r < R2)])] - R2) < 0.8
    # never deeper than -B, and weaker than the wells mid-shell
    assert e.min() >= -B * 2 - 1e-6
    assert e[np.argmin(abs(r - mid))] > e.min() + 0.5 * B * (0.5 if mode == "logistic_shell" else 1)


def test_unknown_mode_raises():
    with pytest.raises(ValueError):
        _energies("nope", [3.0])

"""Initial-structure generators and small coordinate/format helpers."""
import numpy as np
import pytest

from multimm import initial_structure_tools as ist
from multimm.enums import InitialStructureType


@pytest.mark.parametrize("fn", [
    ist.helix_structure, ist.spiral_structure, ist.sphere_surface_structure,
    ist.confined_random_walk, ist.trefoil_knot_structure, ist.random_walk_structure,
])
def test_generators_shape(fn):
    V = fn(50)
    assert V.shape == (50, 3) and np.isfinite(V).all()


@pytest.mark.parametrize("mode", list(InitialStructureType))
def test_compute_init_struct_all_modes(mode):
    V = ist.compute_init_struct(40, mode=mode)
    assert V.shape[1] == 3 and len(V) == 40


def test_polymer_circle():
    V = ist.polymer_circle(20)
    assert np.asarray(V).shape[0] == 20


def test_random_walk_and_saw():
    assert ist.random_walk(30).shape == (30, 3)
    saw = ist.self_avoiding_random_walk(15)
    d = np.linalg.norm(saw[:, None] - saw[None], axis=-1)
    assert d[np.triu_indices(15, 1)].min() > 0.5


def test_generate_hilbert_curve():
    assert ist.generate_hilbert_curve(64, p=3).shape[1] == 3


def test_find_element_indexes():
    assert ist.find_element_indexes(np.array([0, 1, 1, 1, 0, 1]), 1) == "1-3, 5"


def test_hybrid36_roundtrip():
    for val in [0, 5, 9999, 10000, 123456]:
        assert ist.hy36decode(5, ist.hy36encode(5, val)) == val


def test_encode_decode_pure():
    digits = "0123456789"
    assert ist.decode_pure({c: i for i, c in enumerate(digits)}, ist.encode_pure(digits, 4711)) == 4711

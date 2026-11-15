"""Loop import: bead conversion and strength-to-distance rule."""
import numpy as np
import pytest

from multimm.utils import import_mns_from_bedpe

LOOPS = "tests/fixtures/ENCFF045MJY_simple.bedpe"
N_BEADS = 1000


@pytest.fixture()
def loops(tmp_path):
    (tmp_path / "metadata").mkdir()
    return import_mns_from_bedpe(
        LOOPS, N_BEADS, coords=[1, 248387328], chrom="chr1", path=str(tmp_path) + "/"
    )


def test_anchors_are_valid_beads(loops):
    ms, ns, ds, _, _ = loops
    assert len(ms) == len(ns) == len(ds) > 0
    assert ms.min() >= 0 and ns.max() < N_BEADS
    assert np.all(ns - ms > 2)  # loops shorter than 3 beads are dropped


def test_no_duplicate_loops(loops):
    ms, ns, _, _, _ = loops
    assert len({(m, n) for m, n in zip(ms, ns)}) == len(ms)


def test_distance_rule(loops):
    _, _, ds, _, _ = loops
    # strength-derived distances lie in [0.1, 0.2] nm; equal strengths give 1.0
    assert (np.all((ds >= 0.1 - 1e-9) & (ds <= 0.2 + 1e-9))) or np.all(ds == 1.0)


def test_downsampling_keeps_a_subset(tmp_path):
    (tmp_path / "metadata").mkdir()
    full = import_mns_from_bedpe(LOOPS, N_BEADS, coords=[1, 248387328], chrom="chr1", path=str(tmp_path) + "/")
    part = import_mns_from_bedpe(
        LOOPS, N_BEADS, coords=[1, 248387328], chrom="chr1", path=str(tmp_path) + "/", down_prob=0.3
    )
    assert 0 < len(part[0]) < len(full[0])

"""Quick model-build checks: data loading, system and force set-up (no minimisation, no MD)."""
import pytest

from multimm import SimulationEngine
from multimm.config import SimulationConfig
from multimm.model import MultiMM

LOOPS = "tests/fixtures/ENCFF045MJY_simple.bedpe"
COMPS = "tests/fixtures/ENCFF784NCU_reformatted.bed"


def _build(tmp_path, **kw):
    cfg = SimulationConfig(LOOPS_PATH=LOOPS, COMPARTMENT_PATH=COMPS, OUT_PATH=str(tmp_path),
                           N_BEADS=300, SIM_RUN_MD=False, SAVE_PLOTS=False,
                           COB_USE_COMPARTMENT_BLOCKS=True, **kw)
    md = MultiMM(cfg)
    md.set_radiuses()
    md.initialize_simulation()
    md.add_forcefield()
    return md


@pytest.mark.parametrize("kw", [
    dict(CHROM="chr1"),
    dict(CHROM="chr6", LOC_START=1, LOC_END=100_000_000),
    dict(CHROM=None),
])
def test_build_system(tmp_path, kw):
    md = _build(tmp_path, **kw)
    assert md.system.getNumParticles() == 300
    assert md.system.getNumForces() >= 3
    assert "compartment" in " ".join(md.force_groups).lower() or len(md.force_groups) >= 3


def test_genome_wide_has_several_chromosomes(tmp_path):
    md = _build(tmp_path, CHROM=None)
    assert len(md.chrom_idxs) > 1


def test_engine_param_validation():
    ok = SimulationEngine.validate_params({"N_BEADS": 100})
    assert ok["N_BEADS"] == 100
    with pytest.raises(Exception):
        SimulationEngine.validate_params({"NOT_A_KEY": 1})

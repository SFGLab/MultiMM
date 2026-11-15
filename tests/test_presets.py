"""MODELLING_LEVEL presets."""
import pytest

from multimm.config import SimulationConfig
from multimm.run import ArgumentChanger
from multimm.utils import chrom_sizes

LOOPS = "tests/fixtures/ENCFF045MJY_simple.bedpe"
COMPS = "tests/fixtures/ENCFF784NCU_reformatted.bed"


def apply(level, **kwargs):
    cfg = SimulationConfig(LOOPS_PATH=LOOPS, MODELLING_LEVEL=level, **kwargs)
    ArgumentChanger(cfg, chrom_sizes).convenient_argument_changer()
    return cfg


@pytest.mark.parametrize(
    "level, n_beads, container, run_md",
    [("GENE", 1000, False, True), ("REGION", 5000, False, True), ("GW", 200000, True, False)],
)
def test_preset_values(level, n_beads, container, run_md):
    cfg = apply(level)
    assert cfg.N_BEADS == n_beads
    assert cfg.SC_USE_SPHERICAL_CONTAINER is container
    assert cfg.SIM_RUN_MD is run_md
    assert not cfg.CF_USE_CENTRAL_FORCE
    assert not cfg.CHB_USE_CHROMOSOMAL_BLOCKS


def test_chromosome_preset_sets_region():
    cfg = apply("CHROM", CHROM="chr1")
    assert cfg.N_BEADS == 20000
    assert (cfg.LOC_START, cfg.LOC_END) == (1, chrom_sizes["chr1"])


def test_compartment_force_follows_compartment_file():
    assert apply("REGION", COMPARTMENT_PATH=COMPS).COB_USE_COMPARTMENT_BLOCKS
    assert not apply("REGION").COB_USE_COMPARTMENT_BLOCKS


def test_genome_wide_lamina_needs_compartments():
    assert apply("GW", COMPARTMENT_PATH=COMPS).IBL_USE_B_LAMINA_INTERACTION
    assert not apply("GW").IBL_USE_B_LAMINA_INTERACTION


def test_no_level_keeps_user_values():
    cfg = apply("", N_BEADS=1234)
    assert cfg.N_BEADS == 1234

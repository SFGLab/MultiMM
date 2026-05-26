import pytest
from pydantic import ValidationError

from multimm.config import SimulationConfig
from multimm.run import args_tests


def make_config(**kwargs) -> SimulationConfig:
    defaults = dict(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        OUT_PATH="/tmp/output",
    )
    defaults.update(kwargs)
    return SimulationConfig(**defaults)


class TestRequiredPaths:
    def test_missing_loops_raises(self):
        with pytest.raises(ValidationError):
            SimulationConfig(LOOPS_PATH=None, OUT_PATH="/tmp/output")

    def test_empty_loops_raises(self):
        with pytest.raises(ValidationError):
            SimulationConfig(LOOPS_PATH="", OUT_PATH="/tmp/output")

    def test_args_tests_with_valid_path_passes(self):
        cfg = make_config(LOOPS_PATH="/some/file.bedpe")
        args_tests(cfg)


class TestCompartmentConflicts:
    def test_compartment_blocks_without_bed_raises(self):
        cfg = make_config(COB_USE_COMPARTMENT_BLOCKS=True)
        with pytest.raises(ValueError, match="COB_USE_COMPARTMENT_BLOCKS"):
            args_tests(cfg)

    def test_subcompartment_blocks_without_bed_raises(self):
        cfg = make_config(SCB_USE_SUBCOMPARTMENT_BLOCKS=True)
        with pytest.raises(ValueError, match="SCB_USE_SUBCOMPARTMENT_BLOCKS"):
            args_tests(cfg)

    def test_lamina_without_compartments_raises(self):
        cfg = make_config(IBL_USE_B_LAMINA_INTERACTION=True)
        with pytest.raises(ValueError):
            args_tests(cfg)

    def test_lamina_without_compartment_force_raises(self):
        cfg = make_config(
            IBL_USE_B_LAMINA_INTERACTION=True,
            COMPARTMENT_PATH="/some/comps.bed",
        )
        with pytest.raises(ValueError):
            args_tests(cfg)

    def test_lamina_with_compartment_force_passes(self):
        cfg = make_config(
            IBL_USE_B_LAMINA_INTERACTION=True,
            COMPARTMENT_PATH="/some/comps.bed",
            COB_USE_COMPARTMENT_BLOCKS=True,
        )
        args_tests(cfg)


class TestNucleosomeConflicts:
    def test_nuc_interpolation_without_atacseq_raises(self):
        cfg = make_config(NUC_DO_INTERPOLATION=True)
        with pytest.raises(ValueError, match="NUC_DO_INTERPOLATION"):
            args_tests(cfg)

    def test_nuc_interpolation_with_atacseq_passes(self):
        cfg = make_config(NUC_DO_INTERPOLATION=True, ATACSEQ_PATH="/some/track.bw")
        args_tests(cfg)


class TestChromConflicts:
    def test_central_force_with_single_chrom_raises(self):
        cfg = make_config(CF_USE_CENTRAL_FORCE=True, CHROM="chr1")
        with pytest.raises(ValueError, match="CF_USE_CENTRAL_FORCE"):
            args_tests(cfg)

    def test_chromosomal_blocks_with_single_chrom_raises(self):
        cfg = make_config(CHB_USE_CHROMOSOMAL_BLOCKS=True, CHROM="chr1")
        with pytest.raises(ValueError, match="CHB_USE_CHROMOSOMAL_BLOCKS"):
            args_tests(cfg)

    def test_central_force_genome_wide_passes(self):
        cfg = make_config(CF_USE_CENTRAL_FORCE=True)
        args_tests(cfg)

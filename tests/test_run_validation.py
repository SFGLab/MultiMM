"""Input checks done by `args_tests` (errors) and advisory warnings."""
import logging

import pytest

from multimm.config import SimulationConfig
from multimm.run import args_tests

LOOPS = "tests/fixtures/ENCFF045MJY_simple.bedpe"
COMPS = "tests/fixtures/ENCFF784NCU_reformatted.bed"


def make_config(**kwargs) -> SimulationConfig:
    params = dict(LOOPS_PATH=LOOPS, OUT_PATH="/tmp/output")
    params.update(kwargs)
    return SimulationConfig(**params)


@pytest.fixture()
def track(tmp_path) -> str:
    """Empty file standing in for a .bw track (only existence is checked)."""
    path = tmp_path / "track.bw"
    path.write_bytes(b"")
    return str(path)


class TestRequiredInputs:
    def test_loops_required_by_default(self):
        # LE_USE_HARMONIC_BOND is on by default, so loops are needed
        with pytest.raises(ValueError, match="LOOPS_PATH"):
            args_tests(SimulationConfig(OUT_PATH="/tmp/output"))

    def test_loops_not_needed_when_loop_force_off(self):
        args_tests(SimulationConfig(OUT_PATH="/tmp/output", LE_USE_HARMONIC_BOND=False))

    def test_missing_file_is_reported(self):
        with pytest.raises(ValueError, match="not found"):
            args_tests(make_config(LOOPS_PATH="/some/file.bedpe"))

    def test_valid_loops_pass(self):
        args_tests(make_config())

    def test_hic_force_requires_path(self):
        with pytest.raises(ValueError, match="HIC_PATH"):
            args_tests(make_config(HIC_USE_FORCE=True))


class TestCompartmentConflicts:
    def test_compartment_blocks_without_bed_raises(self):
        with pytest.raises(ValueError, match="COB_USE_COMPARTMENT_BLOCKS"):
            args_tests(make_config(COB_USE_COMPARTMENT_BLOCKS=True))

    def test_subcompartment_blocks_without_bed_raises(self):
        with pytest.raises(ValueError, match="SCB_USE_SUBCOMPARTMENT_BLOCKS"):
            args_tests(make_config(SCB_USE_SUBCOMPARTMENT_BLOCKS=True))

    def test_lamina_without_compartments_raises(self):
        with pytest.raises(ValueError, match="compartment"):
            args_tests(make_config(IBL_USE_B_LAMINA_INTERACTION=True))

    def test_lamina_without_compartment_force_raises(self):
        cfg = make_config(IBL_USE_B_LAMINA_INTERACTION=True, COMPARTMENT_PATH=COMPS)
        with pytest.raises(ValueError, match="compartment-based forces"):
            args_tests(cfg)

    @pytest.mark.parametrize("flag", ["COB_USE_COMPARTMENT_BLOCKS", "SCB_USE_SUBCOMPARTMENT_BLOCKS"])
    def test_lamina_with_compartment_force_passes(self, flag):
        args_tests(make_config(IBL_USE_B_LAMINA_INTERACTION=True, COMPARTMENT_PATH=COMPS, **{flag: True}))


class TestNucleosomeConflicts:
    def test_interpolation_without_atacseq_raises(self):
        with pytest.raises(ValueError, match="NUC_DO_INTERPOLATION"):
            args_tests(make_config(NUC_DO_INTERPOLATION=True))

    def test_interpolation_with_atacseq_passes(self, track):
        args_tests(make_config(NUC_DO_INTERPOLATION=True, ATACSEQ_PATH=track))


class TestAdvisoryWarnings:
    """Questionable but valid setups only warn (they must not raise)."""

    @pytest.mark.parametrize(
        "kwargs, text",
        [
            (dict(CF_USE_CENTRAL_FORCE=True, CHROM="chr1"), "Central force"),
            (dict(CHB_USE_CHROMOSOMAL_BLOCKS=True, CHROM="chr1"), "Chromosomal block"),
            (dict(SC_USE_SPHERICAL_CONTAINER=True, CHROM="chr1"), "spherical container"),
            (dict(EV_USE_EXCLUDED_VOLUME=False), "backbone"),
        ],
    )
    def test_warns_but_passes(self, caplog, kwargs, text):
        with caplog.at_level(logging.WARNING, logger="multimm.run"):
            args_tests(make_config(**kwargs))
        assert text.lower() in caplog.text.lower()

    def test_genome_wide_central_force_is_silent(self, caplog):
        with caplog.at_level(logging.WARNING, logger="multimm.run"):
            args_tests(make_config(CF_USE_CENTRAL_FORCE=True))
        assert "Central force" not in caplog.text


def test_file_extensions_are_case_insensitive(tmp_path):
    f = tmp_path / "comps.bigWig"
    f.write_text("x")
    cfg = make_config(COMPARTMENT_PATH=str(f))
    args_tests(cfg)  # .bigWig accepted
    # the path is found even if its letter case differs from the file on disk
    cfg = make_config(COMPARTMENT_PATH=str(tmp_path / "COMPS.BIGWIG"))
    args_tests(cfg)
    assert cfg.COMPARTMENT_PATH.endswith("comps.bigWig")


def test_wrong_extension_and_missing_file_are_told_apart(tmp_path):
    bad = tmp_path / "comps.txt"
    bad.write_text("x")
    with pytest.raises(ValueError, match="unsupported extension"):
        args_tests(make_config(COMPARTMENT_PATH=str(bad)))
    with pytest.raises(ValueError, match="not found"):
        args_tests(make_config(COMPARTMENT_PATH=str(tmp_path / "nope.bed")))


def test_module_runs_without_runpy_warning():
    import subprocess, sys
    r = subprocess.run([sys.executable, "-W", "error", "-m", "multimm.run", "--help"],
                       capture_output=True, text=True)
    assert "RuntimeWarning" not in r.stderr


def test_ensemble_needs_a_size():
    with pytest.raises(ValueError, match="N_ENSEMBLE"):
        args_tests(make_config(GENERATE_ENSEMBLE=True))
    args_tests(make_config(GENERATE_ENSEMBLE=True, N_ENSEMBLE=3))

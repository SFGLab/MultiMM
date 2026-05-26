import os
from multimm.config import SimulationConfig
from multimm.model import MultiMM


def test_simulation_chrom1(tmp_path):
    out_dir = tmp_path / "sim_chrom1"
    config = SimulationConfig(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        COMPARTMENT_PATH="tests/fixtures/ENCFF784NCU_reformatted.bed",
        OUT_PATH=str(out_dir),
        N_BEADS=1000,
        CHROM="chr1",
        LOC_START=1,
        LOC_END=248387328,
        SIM_RUN_MD=True,
        SIM_N_STEPS=10,
        SAVE_PLOTS=False,
        COB_USE_COMPARTMENT_BLOCKS=True,
    )
    md = MultiMM(config)
    md.run()
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))


def test_simulation_chrom1_no_coords(tmp_path):
    out_dir = tmp_path / "sim_chrom1_no_coords"
    config = SimulationConfig(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        COMPARTMENT_PATH="tests/fixtures/ENCFF784NCU_reformatted.bed",
        OUT_PATH=str(out_dir),
        N_BEADS=1000,
        CHROM="chr1",
        SIM_RUN_MD=True,
        SIM_N_STEPS=10,
        SAVE_PLOTS=False,
        COB_USE_COMPARTMENT_BLOCKS=True,
    )
    md = MultiMM(config)
    md.run()
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))



def test_simulation_genome_wide(tmp_path):
    out_dir = tmp_path / "sim_gw"
    config = SimulationConfig(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        COMPARTMENT_PATH="tests/fixtures/ENCFF784NCU_reformatted.bed",
        OUT_PATH=str(out_dir),
        N_BEADS=1000,
        CHROM=None,
        SIM_RUN_MD=True,
        SIM_N_STEPS=10,
        SAVE_PLOTS=False,
        COB_USE_COMPARTMENT_BLOCKS=True,
    )
    md = MultiMM(config)
    md.run()
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))


def test_engine_in_process_run(tmp_path):
    out_dir = tmp_path / "bridge_test"
    params = {
        "LOOPS_PATH": "tests/fixtures/ENCFF045MJY_simple.bedpe",
        "COMPARTMENT_PATH": "tests/fixtures/ENCFF784NCU_reformatted.bed",
        "OUT_PATH": str(out_dir),
        "N_BEADS": 1000,
        "CHROM": "1",
        "LOC_START": 1,
        "LOC_END": 248387328,
        "SIM_RUN_MD": True,
        "SIM_N_STEPS": 10,
        "SAVE_PLOTS": False,
        "COB_USE_COMPARTMENT_BLOCKS": True,
    }
    from multimm import SimulationEngine
    config_path = SimulationEngine.run_in_process(params)
    assert os.path.exists(config_path)
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))
    assert os.path.exists(os.path.join(str(out_dir), "metadata", "output.log"))


def test_subprocess_run(tmp_path):
    out_dir = tmp_path / "subproc_test"
    params = {
        "LOOPS_PATH": "tests/fixtures/ENCFF045MJY_simple.bedpe",
        "COMPARTMENT_PATH": "tests/fixtures/ENCFF784NCU_reformatted.bed",
        "OUT_PATH": str(out_dir),
        "N_BEADS": 1000,
        "CHROM": "1",
        "LOC_START": 1,
        "LOC_END": 248387328,
        "SIM_RUN_MD": True,
        "SIM_N_STEPS": 10,
        "SAVE_PLOTS": False,
        "COB_USE_COMPARTMENT_BLOCKS": True,
    }
    from multimm import SimulationEngine
    config_path = SimulationEngine.run_subprocess(params)
    assert os.path.exists(config_path)
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))
    assert os.path.exists(os.path.join(str(out_dir), "metadata", "output.log"))



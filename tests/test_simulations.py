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


def test_simulation_chrom6_names(tmp_path):
    out_dir = tmp_path / "sim_chrom6_names"
    config = SimulationConfig(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        COMPARTMENT_PATH="tests/fixtures/ENCFF784NCU_reformatted.bed",
        OUT_PATH=str(out_dir),
        N_BEADS=1000,
        CHROM="chr6",
        SIM_RUN_MD=True,
        SIM_N_STEPS=10,
        SAVE_PLOTS=True,
        COB_USE_COMPARTMENT_BLOCKS=True,
    )
    md = MultiMM(config)
    md.run()
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))
    assert os.path.exists(os.path.join(str(out_dir), "model", "chromosomes", "MultiMM_minimized_chr6.cif"))
    assert not os.path.exists(os.path.join(str(out_dir), "model", "chromosomes", "MultiMM_minimized_chr1.cif"))
    assert os.path.exists(os.path.join(str(out_dir), "plots", "chromosomes", "chr6_minimized_structure.png"))
    assert not os.path.exists(os.path.join(str(out_dir), "plots", "chromosomes", "chr1_minimized_structure.png"))




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


def test_simulation_region_no_compartments_with_plots(tmp_path):
    out_dir = tmp_path / "sim_region_no_comps"
    config = SimulationConfig(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        COMPARTMENT_PATH=None,
        OUT_PATH=str(out_dir),
        N_BEADS=100,
        CHROM="chr1",
        LOC_START=16000000,
        LOC_END=17950000,
        SIM_RUN_MD=True,
        SIM_N_STEPS=5,
        SAVE_PLOTS=True,
        COB_USE_COMPARTMENT_BLOCKS=False,
    )
    md = MultiMM(config)
    md.run()
    assert os.path.exists(os.path.join(str(out_dir), "model", "MultiMM_minimized.cif"))




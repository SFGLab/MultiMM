import textwrap

import pytest

from multimm.config import SimulationConfig


@pytest.fixture()
def default_config() -> SimulationConfig:
    return SimulationConfig()


@pytest.fixture()
def minimal_config() -> SimulationConfig:
    return SimulationConfig(
        LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
        OUT_PATH="/tmp/multimm_output",
    )


@pytest.fixture()
def sample_ini(tmp_path) -> str:
    ini_content = textwrap.dedent("""\
        [Main]
        PLATFORM = CPU
        N_BEADS = 1000
        LOOPS_PATH = tests/fixtures/ENCFF045MJY_simple.bedpe
        OUT_PATH = /tmp/output
        SIM_RUN_MD = False
        SIM_N_STEPS = 500
    """)
    ini_file = tmp_path / "test_config.ini"
    ini_file.write_text(ini_content)
    return str(ini_file)

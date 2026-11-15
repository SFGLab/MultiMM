"""Configuration: defaults, parsing and key checking."""
import configparser

import pytest
from openmm.unit import kelvin, nanometer
from pydantic import ValidationError

from multimm.config import SimulationConfig
from multimm.run import _check_known_config_fields, my_config_parser

# Defaults that the documentation quotes (docs/model/*.md)
DOCUMENTED_DEFAULTS = {
    "N_BEADS": 50000,
    "POL_USE_HARMONIC_BOND": True,
    "POL_USE_HARMONIC_ANGLE": True,
    "EV_USE_EXCLUDED_VOLUME": True,
    "EV_FORCE_TYPE": "powerlaw",
    "EV_EPSILON": 100.0,
    "EV_POWER": 6.0,
    "EV_R_SMALL": 0.05,
    "LE_USE_HARMONIC_BOND": True,
    "LE_LOOP_FORCE_TYPE": "harmonic",
    "LE_FIXED_DISTANCES": False,
    "COB_USE_COMPARTMENT_BLOCKS": False,
    "COB_FORCE_TYPE": "gaussian",
    "COB_EA": 1.0,
    "COB_EB": 2.0,
    "SCB_FORCE_TYPE": "gaussian",
    "SCB_EA1": 1.0,
    "SCB_EA2": 1.33,
    "SCB_EB1": 1.66,
    "SCB_EB2": 2.0,
    "CHB_FORCE_TYPE": "polynomial",
    "CHB_KC": 0.3,
    "CHB_DE": 1e-4,
    "SC_SCALE": 1000.0,
    "IBL_SCALE": 400.0,
    "BLAMINA_FORCE_TYPE": "sin",
    "CF_STRENGTH": 20.0,
    "CENTRAL_FORCE_TYPE": "harmonic",
    "HIC_USE_FORCE": False,
    "HIC_K_SCALE": 40.0,
    "HIC_BOLTZMANN_ALPHA": 4.0,
    "HIC_BOLTZMANN_KERNEL": "exponential",
    "HIC_NORMALIZATION": "KR",
    "HIC_AUTO_SCALE": True,
    "HIC_BLOCK_COPOLYMER": False,
    "SIM_RUN_MD": False,
    "SIM_INTEGRATOR_TYPE": "langevin",
    "SIM_FRICTION_COEFF": 0.5,
    "NUC_DO_INTERPOLATION": False,
}


@pytest.mark.parametrize("name, value", DOCUMENTED_DEFAULTS.items())
def test_documented_defaults(default_config, name, value):
    assert getattr(default_config, name) == pytest.approx(value)


def test_quantity_defaults(default_config):
    assert default_config.POL_HARMONIC_BOND_R0.value_in_unit(nanometer) == pytest.approx(0.1)
    assert default_config.LE_HARMONIC_BOND_R0.value_in_unit(nanometer) == pytest.approx(0.1)
    assert default_config.SIM_TEMPERATURE.value_in_unit(kelvin) == pytest.approx(310.0)


def test_quantity_from_text():
    cfg = SimulationConfig(POL_HARMONIC_BOND_R0="0.2 nanometer")
    assert cfg.POL_HARMONIC_BOND_R0.value_in_unit(nanometer) == pytest.approx(0.2)


@pytest.mark.parametrize("bad", ["0.2", 0.2, "nanometer"])
def test_quantity_requires_unit(bad):
    with pytest.raises(ValidationError):
        SimulationConfig(POL_HARMONIC_BOND_R0=bad)


@pytest.mark.parametrize("text, expected", [("yes", True), ("1", True), ("False", False), ("no", False)])
def test_boolean_spellings(text, expected):
    assert SimulationConfig(SIM_RUN_MD=text).SIM_RUN_MD is expected


@pytest.mark.parametrize("given, expected", [("1", "chr1"), ("chr6", "chr6"), ("", None), ("none", None)])
def test_chromosome_normalised(given, expected):
    assert SimulationConfig(CHROM=given).CHROM == expected


def test_unknown_key_rejected():
    with pytest.raises(ValidationError):
        SimulationConfig(NOT_A_REAL_KEY=1)


def test_known_key_check_suggests_close_name():
    parser = configparser.ConfigParser()
    parser.read_string("[Main]\nHIC_BOLTZMAN_ALPHA = 4\n")
    with pytest.raises(ValueError, match="HIC_BOLTZMANN_ALPHA"):
        _check_known_config_fields(my_config_parser(parser), "x.ini")


def test_known_keys_accepted():
    parser = configparser.ConfigParser()
    parser.read_string("[Main]\nN_BEADS = 100\nSIM_RUN_MD = True\n")
    _check_known_config_fields(my_config_parser(parser), "x.ini")

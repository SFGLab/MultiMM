"""Documented defects. Each test is expected to fail (xfail, strict) until
the defect is fixed; when it is fixed pytest reports XPASS and the marker
should be removed."""
import inspect
import re

import pytest

from multimm import model
from multimm.nucleosome_interpolation import NucleosomeInterpolation


@pytest.mark.xfail(strict=True, reason="model.nuc_interpolation passes 6 positional args to a 5-arg __init__")
def test_nucleosome_call_matches_signature():
    src = inspect.getsource(model.MultiMM.nuc_interpolation)
    n_passed = len(re.search(r"NucleosomeInterpolation\((.*?)\n\s*\)", src, re.S).group(1).split(","))
    n_params = len(inspect.signature(NucleosomeInterpolation.__init__).parameters) - 1  # self
    assert n_passed <= n_params + 1  # trailing comma gives one empty item


@pytest.mark.xfail(strict=True, reason="yukawa compartment energy uses s1 twice, so it does not segregate A from B")
def test_yukawa_compartment_uses_both_beads():
    src = inspect.getsource(model.MultiMM.add_compartment_blocks)
    yukawa = src.split('mode == "yukawa"')[1].split('mode == "theta"')[0]
    assert "delta(s2-1)" in yukawa

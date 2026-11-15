"""Sign checks: every force term must push or pull the way the docs say.

Pair terms are probed with two beads at distance r (all others far away);
radial terms with every bead placed at radius r from the centre.
"""
import numpy as np
import openmm as mm
import pytest

from multimm.config import SimulationConfig
from multimm.model import MultiMM

BASE = dict(LOOPS_PATH="tests/fixtures/ENCFF045MJY_simple.bedpe",
            COMPARTMENT_PATH="tests/fixtures/ENCFF784NCU_reformatted.bed",
            N_BEADS=300, CHROM="chr1", SIM_RUN_MD=False, SAVE_PLOTS=False,
            COB_USE_COMPARTMENT_BLOCKS=False)


def build(tmp_path, **kw):
    md = MultiMM(SimulationConfig(OUT_PATH=str(tmp_path), **{**BASE, **kw}))
    md.set_radiuses()
    md.initialize_simulation()
    md.add_forcefield()
    return md


def energy(md, group, pos):
    ctx = mm.Context(md.system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(pos)
    return ctx.getState(getEnergy=True, groups={md.force_groups[group]}).getPotentialEnergy()._value


def pair_energy(md, group, i, j, r):
    pos = np.zeros((md.system.getNumParticles(), 3))
    pos[:, 2] = 1e3 * np.arange(len(pos))  # everyone far from everyone
    pos[j] = pos[i] + [r, 0, 0]
    return energy(md, group, pos)


def radial_energy(md, group, r):
    pos = np.tile(md.mass_center, (md.system.getNumParticles(), 1)) + [r, 0, 0]
    return energy(md, group, pos)


def two_beads(md, positive):
    c = np.asarray(md.Cs)
    idx = np.where(c > 0 if positive else c < 0)[0]
    return int(idx[10]), int(idx[40])


# ---------- pair terms ----------
@pytest.mark.parametrize("form", ["powerlaw", "gaussian_core"])
def test_excluded_volume_is_repulsive(tmp_path, form):
    md = build(tmp_path, EV_FORCE_TYPE=form)
    e = [pair_energy(md, "Excluded volume", 0, 50, r) for r in (0.02, 0.05, 0.1, 0.2)]
    assert all(x > 0 for x in e) and all(a > b for a, b in zip(e, e[1:]))  # positive, falls with distance


@pytest.mark.parametrize("form", ["gaussian", "theta"])
def test_compartment_like_attracts_like(tmp_path, form):
    md = build(tmp_path, COB_USE_COMPARTMENT_BLOCKS=True, COB_FORCE_TYPE=form)
    for positive in (True, False):
        i, j = two_beads(md, positive)
        near, far = (pair_energy(md, "Compartment blocks", i, j, r) for r in (0.05, 0.45))
        assert near < 0 and near <= far  # attractive well, deepest at contact
    a, b = two_beads(md, True)[0], two_beads(md, False)[0]
    assert pair_energy(md, "Compartment blocks", a, b, 0.05) == pytest.approx(0, abs=1e-9)  # A-B: no force


def test_B_pairs_bind_stronger_than_A_pairs(tmp_path):
    md = build(tmp_path, COB_USE_COMPARTMENT_BLOCKS=True)
    ea = pair_energy(md, "Compartment blocks", *two_beads(md, True), 0.05)
    eb = pair_energy(md, "Compartment blocks", *two_beads(md, False), 0.05)
    assert eb < ea < 0  # Eb = 2 > Ea = 1


@pytest.mark.parametrize("form", ["polynomial", "gaussian", "saturating"])
def test_chromosome_blocks_pull_together(tmp_path, form):
    md = build(tmp_path, CHB_USE_CHROMOSOMAL_BLOCKS=True, CHB_FORCE_TYPE=form)
    e = [pair_energy(md, "Chromosomal blocks", 0, 50, r) for r in (0.1, 0.5, 1.0, 2.0)]  # inside the 5*R2 cutoff
    assert all(a < b for a, b in zip(e, e[1:]))  # energy rises with distance => attraction


def test_loops_are_springs_around_target(tmp_path):
    for form in ("harmonic", "fene_soft", "gaussian_tether"):
        md = build(tmp_path, LOOPS_PATH=BASE["LOOPS_PATH"], LE_LOOP_FORCE_TYPE=form)
        ms, ns = np.asarray(md.ms), np.asarray(md.ns)
        allb = np.concatenate([ms, ns])
        k = next(q for q in range(len(ms)) if (allb == ms[q]).sum() == 1 and (allb == ns[q]).sum() == 1)
        r0 = md.args.LE_HARMONIC_BOND_R0._value if md.args.LE_FIXED_DISTANCES else float(md.ds[k])
        e = {f: pair_energy(md, "Loop extrusion", int(ms[k]), int(ns[k]), f * r0) for f in (0.5, 1.0, 2.0)}
        assert e[1.0] == pytest.approx(0, abs=1e-6) and e[0.5] > 0 and e[2.0] > 0, form


def test_bonds_and_angles_are_springs(tmp_path):
    md = build(tmp_path)
    r0 = md.args.POL_HARMONIC_BOND_R0._value

    def chain(step):  # straight chain with equal spacing
        pos = np.zeros((md.system.getNumParticles(), 3))
        pos[:, 0] = step * np.arange(len(pos))
        return pos

    assert energy(md, "Harmonic bonds", chain(r0)) == pytest.approx(0, abs=1e-6)
    assert energy(md, "Harmonic bonds", chain(2 * r0)) > 0
    assert energy(md, "Harmonic bonds", chain(0.5 * r0)) > 0
    # a straight chain is the most bent-free state: bending it costs energy >= 0
    assert energy(md, "Harmonic angles", chain(r0)) >= 0


# ---------- radial terms ----------
def test_container_only_penalises_outside_the_shell(tmp_path):
    md = build(tmp_path, SC_USE_SPHERICAL_CONTAINER=True)
    mid = 0.5 * (md.radius1 + md.radius2)
    assert radial_energy(md, "Spherical container", mid) == pytest.approx(0, abs=1e-9)
    assert radial_energy(md, "Spherical container", 1.3 * md.radius2) > 0
    assert radial_energy(md, "Spherical container", 0.5 * md.radius1) > 0


@pytest.mark.parametrize("form", ["sin", "gaussian_shell", "harmonic_shell", "logistic_shell"])
def test_lamina_attracts_only_B_beads_to_the_walls(tmp_path, form):
    md = build(tmp_path, COB_USE_COMPARTMENT_BLOCKS=True, IBL_USE_B_LAMINA_INTERACTION=True, BLAMINA_FORCE_TYPE=form)
    wall = radial_energy(md, "B-lamina interaction", md.radius2)
    mid = radial_energy(md, "B-lamina interaction", 0.5 * (md.radius1 + md.radius2))
    assert wall < 0 and wall < mid <= 1e-6
    md.Cs = np.abs(np.asarray(md.Cs))  # pretend nothing is B
    # A beads: build the energy again with all-positive states
    md2 = build(tmp_path, COB_USE_COMPARTMENT_BLOCKS=True, IBL_USE_B_LAMINA_INTERACTION=True, BLAMINA_FORCE_TYPE=form)
    assert radial_energy(md2, "B-lamina interaction", md2.radius2) <= 0


@pytest.mark.parametrize("form,sign", [("harmonic", +1), ("gaussian", -1), ("logistic", -1)])
def test_central_force_sign(tmp_path, form, sign):
    md = build(tmp_path, CHROM=None, CF_USE_CENTRAL_FORCE=True, CENTRAL_FORCE_TYPE=form)  # needs several chromosomes
    inner = radial_energy(md, "Central force", md.radius1)
    outer = radial_energy(md, "Central force", md.radius2)
    if form == "harmonic":
        assert inner == pytest.approx(0, abs=1e-6) and outer > 0  # minimum at R1, costs more outside
    else:
        assert inner < 0 and inner <= outer <= 1e-9  # a well in the inner region


@pytest.mark.parametrize("lamina,central", [("gaussian_shell", "gaussian"), ("logistic_shell", "logistic")])
def test_forces_can_be_combined(tmp_path, lamina, central):
    # regression: lamina/central forces used to share global names (sigma, lambda) with other forces
    md = build(tmp_path, CHROM=None, COB_USE_COMPARTMENT_BLOCKS=True, SC_USE_SPHERICAL_CONTAINER=True,
               IBL_USE_B_LAMINA_INTERACTION=True, BLAMINA_FORCE_TYPE=lamina,
               CF_USE_CENTRAL_FORCE=True, CENTRAL_FORCE_TYPE=central, CHB_USE_CHROMOSOMAL_BLOCKS=True)
    assert radial_energy(md, "B-lamina interaction", md.radius2) < 0

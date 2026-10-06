"""
hic_force.py — builds the OpenMM Hi-C contact force.

Boltzmann-inversion PMF: each pair's contact strength c_ij is converted to
a target distance r_target via one of three P(r) kernels (see
VALID_BOLTZMANN_KERNELS), then restrained there:

    U_ij(r) = 0.5 * k_scale * c_ij * (r - r_target)^2        (real evidence)
    U_ij(r) = 0.5 * k_scale * c_ij * max(0, rc - r)^2         (background pair)

Pairs whose kernel saturates at rc carry no distance information beyond
"not enriched" — restraining the (often >50%) majority of such pairs to one
shared exact distance with a two-sided well is what spreads beads onto a
spherical shell (same mechanism as the Thomson problem). So those pairs
automatically get a one-sided floor only (never pulled together, just kept
from overlapping); pairs with real sub-rc evidence keep the full two-sided
well. This split is derived per-pair straight from the data (whether its own
r_target hit the rc cap) — there is no separate tunable for it.
Matrix preprocessing (clean/resize/balance/OE/denoise) lives in
read_hic.preprocess_hic_matrix — this module only builds the force.

Usage (model.py):
    from hic_force import build_hic_force
    force = build_hic_force(H_raw, N_beads=self.N, rc=self.hic_rc,
                             alpha=4.0, k_scale=10.0)
    self.system.addForce(force)
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

try:
    import openmm as mm
except ImportError:                     # legacy simtk namespace
    from simtk import openmm as mm

log = logging.getLogger(__name__)

from .logger import log_table  # noqa: E402 — after log is set up
from .read_hic import preprocess_hic_matrix  # noqa: E402


# ═════════════════════════════════════════════════════════════════════════════
# Distance <-> contact-probability inversion
# ═════════════════════════════════════════════════════════════════════════════

# Three P(r) kernels for c_ij -> r_target; alpha is the steepness knob for
# all three. Default "exponential" is the classic Boltzmann distribution.
VALID_BOLTZMANN_KERNELS = ("power_law", "exponential", "sigmoid")


def _boltzmann_r_target(
    c_vals: np.ndarray,
    rc: float,
    alpha: float,
    r_min: float,
    kernel: str = "exponential",
) -> np.ndarray:
    """c_ij -> target distance, clipped to [r_min, rc].

    exponential : r_target = r_min - lam*ln(c),  lam = rc/alpha
    power_law   : r_target = r_min * c^(-1/alpha)
    sigmoid     : r_target = r0 - (1/k)*ln(c/(1-c)),  r0=(r_min+rc)/2, k=alpha/rc

    Each is the exact inverse of its P(r) in `get_boltzmann_p_func`.
    """
    if kernel not in VALID_BOLTZMANN_KERNELS:
        raise ValueError(f"Unknown HIC_BOLTZMANN_KERNEL {kernel!r}; choose from {VALID_BOLTZMANN_KERNELS}")

    if kernel == "power_law":
        c_safe = np.clip(c_vals, 1e-300, 1.0)
        with np.errstate(divide="ignore", over="ignore"):
            r_target = r_min * np.power(c_safe, -1.0 / alpha)
    elif kernel == "exponential":
        lam = rc / alpha
        c_safe = np.clip(c_vals, 1e-300, 1.0)
        r_target = r_min - lam * np.log(c_safe)
    else:  # sigmoid
        r0 = 0.5 * (r_min + rc)
        k = alpha / rc
        c_safe = np.clip(c_vals, 1e-12, 1.0 - 1e-12)
        r_target = r0 - (1.0 / k) * np.log(c_safe / (1.0 - c_safe))

    return np.clip(r_target, r_min, rc)


def get_boltzmann_p_func(
    rc     : float,
    alpha  : float           = None,
    r_min  : Optional[float] = None,
    kernel : str             = "exponential",
):
    """P(r) for the given kernel — the exact inverse of `_boltzmann_r_target`
    (P(r_target(c)) == c). Lets validation.py/plots.py score a structure
    against the same law the force targets, instead of a generic proxy.

    exponential : P(r) = clip(exp(-(r-r_min)/lam), 0, 1), lam = rc/alpha
    power_law   : P(r) = clip((r_min/r)^alpha, 0, 1)
    sigmoid     : P(r) = 1 / (1 + exp(k*(r-r0))), r0=(r_min+rc)/2, k=alpha/rc

    r_min defaults to 0.1*rc. Returns p_func(r) -> P, vectorised.
    """
    if kernel not in VALID_BOLTZMANN_KERNELS:
        raise ValueError(f"Unknown HIC_BOLTZMANN_KERNEL {kernel!r}; choose from {VALID_BOLTZMANN_KERNELS}")
    if alpha is None:
        alpha = _DEFAULT_BOLTZMANN_ALPHA
    if r_min is None or r_min <= 0:
        r_min = 0.1 * rc

    if kernel == "power_law":
        def p_func(r):
            with np.errstate(divide="ignore", over="ignore"):
                p = np.power(r_min / np.maximum(r, 1e-12), alpha)
            return np.clip(p, 0.0, 1.0)
    elif kernel == "exponential":
        lam = rc / alpha
        def p_func(r):
            with np.errstate(over="ignore"):
                p = np.exp(-(np.asarray(r, dtype=np.float64) - r_min) / lam)
            return np.clip(p, 0.0, 1.0)
    else:  # sigmoid
        r0 = 0.5 * (r_min + rc)
        k = alpha / rc
        def p_func(r):
            with np.errstate(over="ignore"):
                p = 1.0 / (1.0 + np.exp(k * (np.asarray(r, dtype=np.float64) - r0)))
            return np.clip(p, 0.0, 1.0)

    return p_func


def auto_contact_scale(
    coords: np.ndarray,
    percentile: float = 10.0,
    n_samples: int = 20000,
    seed: int = 0,
) -> float:
    """Contact length scale from a structure's own pairwise-distance
    distribution (the `percentile`-th of sampled distances), instead of the
    force's microscopic rc/sigma, which underflows for distant pairs.
    """
    coords = np.asarray(coords)
    N = coords.shape[0]
    if N < 2:
        return 1.0
    rng = np.random.default_rng(seed)
    draw = max(n_samples * 2, 64)
    i_idx = rng.integers(0, N, size=draw)
    j_idx = rng.integers(0, N, size=draw)
    keep = i_idx != j_idx
    i_idx, j_idx = i_idx[keep], j_idx[keep]
    n_use = min(n_samples, i_idx.size)
    if n_use == 0:
        return 1.0
    i_idx, j_idx = i_idx[:n_use], j_idx[:n_use]
    d = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)
    d = d[np.isfinite(d) & (d > 0)]
    if d.size == 0:
        return 1.0
    return float(np.percentile(d, percentile))


class HiCNoiseController:
    """Re-noises each bond's c_ij around its original value (r_lo/r_hi stay
    fixed) so sampling explores nearby configurations instead of settling
    on one attractor. Call `resample()` periodically with noise std-dev.
    """

    def __init__(self, force, rows, cols, base_c_vals, extra_vals=None, seed=0):
        self.force = force
        self.rows = np.asarray(rows)
        self.cols = np.asarray(cols)
        self.base_c_vals = np.asarray(base_c_vals, dtype=np.float64)
        # extra_vals: (n_bonds, n_extra) fixed per-bond params appended after c_ij
        self.extra_vals = None if extra_vals is None else np.asarray(extra_vals, dtype=np.float64)
        self.rng = np.random.default_rng(seed)

    def resample(self, context, intensity: float) -> None:
        if intensity <= 0 or self.base_c_vals.size == 0:
            return
        noisy = np.clip(
            self.base_c_vals + self.rng.normal(0.0, intensity, size=self.base_c_vals.shape),
            0.0, 1.0,
        )
        set_params = self.force.setBondParameters
        if self.extra_vals is not None:
            for idx in range(noisy.shape[0]):
                set_params(
                    idx, int(self.rows[idx]), int(self.cols[idx]),
                    [float(noisy[idx]), *self.extra_vals[idx].tolist()],
                )
        else:
            for idx in range(noisy.shape[0]):
                set_params(idx, int(self.rows[idx]), int(self.cols[idx]), [float(noisy[idx])])
        self.force.updateParametersInContext(context)


_DEFAULT_BOLTZMANN_ALPHA = 4.0   # standard Hi-C scaling exponent, c ~ r^-alpha


def build_boltzmann_force(
    C            : np.ndarray,
    rc       : float,
    alpha        : float           = _DEFAULT_BOLTZMANN_ALPHA,
    k_scale      : float           = 1.0,
    force_group  : int             = 1,
    return_controller: bool        = False,
    noise_seed   : int             = 0,
    r_min        : Optional[float] = None,
    kernel       : str             = "exponential",
):
    """Sparse CustomBondForce via Boltzmann-inversion: c_ij -> r_target per
    kernel, then restrained according to whether that target is real
    distance evidence or just hit the kernel's rc cap:

        U_ij(r) = 0.5 * k_scale * c_ij * (r - r_target)^2   (r_target < rc)
        U_ij(r) = 0.5 * k_scale * c_ij * max(0, rc - r)^2    (r_target == rc)

    A pair whose kernel saturates at rc (c_ij too weak/background to pin
    down a real distance) only tells you "not enriched" — not "exactly rc
    away" — so it gets a one-sided floor only: pushed apart if closer than
    rc, free to be anywhere farther. A pair with real sub-rc evidence keeps
    the full two-sided well (pulls if farther, pushes if closer).

    This matters because, with the plain two-sided well applied to every
    pair, weak/background contacts are usually the majority and nearly all
    saturate at the same r_target=rc: restraining that many pairs to one
    shared *exact* distance is only satisfiable in 3-D by spreading beads
    over a spherical shell (same mechanism as the Thomson problem). Exempting
    those pairs from the pull-together side removes that false precision
    automatically — it's derived from each pair's own data, not a tunable.

    Parameters
    ----------
    C           : (N_beads, N_beads) preprocessed contact matrix.
    rc          : upper clip bound [nm] for target distances; also the
                  fallback scale for r_min if not given.
    alpha       : Hi-C scaling exponent, c ~ r^-alpha. Typical 3-4.
    k_scale     : harmonic stiffness [kJ/mol].
    force_group : OpenMM force-group index.
    return_controller : also return a HiCNoiseController if True.
    noise_seed  : RNG seed for the controller.
    r_min       : excluded-volume floor [nm]; falls back to 0.1*rc.
    kernel      : see VALID_BOLTZMANN_KERNELS.

    Returns
    -------
    force, or (force, HiCNoiseController) if return_controller=True.
    """
    if kernel not in VALID_BOLTZMANN_KERNELS:
        raise ValueError(f"Unknown HIC_BOLTZMANN_KERNEL {kernel!r}; choose from {VALID_BOLTZMANN_KERNELS}")

    if r_min is None or r_min <= 0:
        r_min = 0.1 * rc
        log.warning(
            "build_boltzmann_force: r_min not given — using placeholder 0.1*rc=%.4f nm",
            r_min,
        )

    N_beads = C.shape[0]
    log.info(
        "build_boltzmann_force: N=%d, kernel=%s, alpha=%.3f, r_min=%.4f nm, rc=%.3f nm",
        N_beads, kernel, alpha, r_min, rc,
    )

    # normalise to [0, 1]
    C_max = C.max()
    if C_max > 0:
        C = C / C_max
    else:
        log.warning("  contact matrix is all-zero — no bonds will be added")

    # nonzero pairs directly — cheaper than triu_indices+mask at scale (C is
    # guaranteed non-negative, so nonzero == ">0")
    rows_all, cols_all = np.nonzero(C)
    mask       = rows_all < cols_all
    rows, cols = rows_all[mask], cols_all[mask]
    c_vals     = C[rows, cols]

    n_bonds = len(rows)
    log.info("  pairs with nonzero contact strength: %d  (%.2f%% of upper triangle)",
             n_bonds, 100.0 * n_bonds / (N_beads * (N_beads - 1) / 2))

    if n_bonds == 0:
        log.warning("  no bonds added — contact matrix is all-zero after normalisation")

    r_target = _boltzmann_r_target(c_vals, rc, alpha, r_min, kernel=kernel)

    # capped == no real distance evidence (kernel saturated at rc); these
    # get a one-sided floor only (r_hi pushed out past any realistic bead
    # separation), everyone else keeps the exact two-sided well (r_lo==r_hi)
    at_cap = r_target >= rc - 1e-9
    n_at_cap = int(np.count_nonzero(at_cap))
    log.info(
        "  target distances: min=%.4f nm  median=%.4f nm  max=%.4f nm  "
        "(%d/%d pairs, %.1f%%, capped at rc=%.3f nm -> one-sided floor, no pull)",
        float(r_target.min()) if n_bonds else float("nan"),
        float(np.median(r_target)) if n_bonds else float("nan"),
        float(r_target.max()) if n_bonds else float("nan"),
        n_at_cap, n_bonds, 100.0 * n_at_cap / max(n_bonds, 1), rc,
    )

    r_lo = r_target
    r_hi = np.where(at_cap, 10.0 * rc, r_target)

    w_mean = float(c_vals.mean()) if n_bonds else 0.0
    bonds_per_bead = 2.0 * n_bonds / N_beads
    f_typical = k_scale * w_mean * float(np.median(r_target)) * 0.1 if n_bonds else 0.0
    log.info(
        "  calibration: k_scale=%.2f kJ/mol, <c_ij>=%.3f, %.1f bonds/bead, "
        "~%.1f kJ/mol/nm restoring force at 10%% off target",
        k_scale, w_mean, bonds_per_bead, f_typical,
    )
    if k_scale > 160.0:
        log.warning(
            "  HIC_K_SCALE=%.2f kJ/mol is high — may freeze MD sampling (20-80 typical)",
            k_scale,
        )

    # flat-bottom: zero inside [r_lo, r_hi], harmonic beyond either edge
    expression = (
        "0.5 * hic_k * c_ij * (d_below^2 + d_above^2);"
        "d_below = max(0, r_lo - r);"
        "d_above = max(0, r - r_hi)"
    )
    force = mm.CustomBondForce(expression)
    force.addGlobalParameter("hic_k", k_scale)
    force.addPerBondParameter("c_ij")
    force.addPerBondParameter("r_lo")
    force.addPerBondParameter("r_hi")

    # tolist() + local addBond binding avoids per-iteration numpy overhead
    rows_l, cols_l = rows.tolist(), cols.tolist()
    c_vals_l, r_lo_l, r_hi_l = c_vals.tolist(), r_lo.tolist(), r_hi.tolist()
    add_bond = force.addBond
    for r_i, c_i, cv, rlo, rhi in zip(rows_l, cols_l, c_vals_l, r_lo_l, r_hi_l):
        add_bond(r_i, c_i, [cv, rlo, rhi])

    force.setForceGroup(force_group)
    log.info("  Hi-C Boltzmann-PMF force ready  (%d bonds)", n_bonds)

    if not return_controller:
        return force

    controller = HiCNoiseController(
        force, rows, cols, c_vals, extra_vals=np.column_stack([r_lo, r_hi]), seed=noise_seed,
    )
    return force, controller


# ═════════════════════════════════════════════════════════════════════════════
# Public entry point
# ═════════════════════════════════════════════════════════════════════════════

def build_hic_force(
    H_raw            : np.ndarray,
    N_beads          : int,
    rc           : float,
    alpha            : float           = _DEFAULT_BOLTZMANN_ALPHA,
    k_scale          : float           = 1.0,
    force_group      : int             = 1,
    already_balanced : bool            = False,
    oe_normalize     : bool            = False,
    return_controller: bool            = False,
    noise_seed       : int             = 0,
    save_path        : Optional[str]   = None,
    chrom            : Optional[str]   = None,
    r_min            : Optional[float] = None,
    kernel           : str             = "exponential",
):
    """Raw Hi-C array -> sparse Hi-C restraint force.

    Pipeline: read_hic.preprocess_hic_matrix() (clean/resize/balance/OE/
    denoise) -> build_boltzmann_force().

    Parameters
    ----------
    H_raw            : (M, M) raw Hi-C counts; resampled if M != N_beads.
    N_beads          : number of simulation beads.
    rc               : characteristic contact distance [nm].
    alpha            : Hi-C scaling exponent, typical 3-4.
    k_scale          : harmonic stiffness [kJ/mol], typical 20-80.
    force_group      : OpenMM force-group index.
    already_balanced : True -> skip diagonal_normalize().
    oe_normalize     : True -> target OE enrichment above background
                       (see read_hic.oe_enrichment_matrix) instead of raw
                       contact frequency.
    return_controller: also return a HiCNoiseController if True.
    noise_seed       : RNG seed for the controller.
    save_path/chrom  : if given, saves a before/after denoising plot.
    r_min            : excluded-volume floor [nm].
    kernel           : see VALID_BOLTZMANN_KERNELS.

    Returns
    -------
    force, or (force, HiCNoiseController) if return_controller=True.

    Example
    -------
    >>> force = build_hic_force(hic_array, N_beads=500, rc=0.15,
    ...                         alpha=4.0, k_scale=10.0)
    >>> system.addForce(force)
    """
    log_table(
        [
            ("N beads",       N_beads),
            ("kernel",        kernel),
            ("rc",        f"{rc:.3f} nm"),
            ("alpha",         alpha),
            ("k_scale",       f"{k_scale:.3f} kJ/mol"),
            ("Force group",   force_group),
            ("OE normalize",  oe_normalize),
            ("Input shape",   str(H_raw.shape)),
            ("r_min (EV floor)", f"{r_min:.4f} nm" if r_min is not None else "not given — falls back to 0.1*rc"),
        ],
        title="Hi-C Force — build",
        log_fn=log.info,
    )

    H = preprocess_hic_matrix(
        H_raw, N_beads, already_balanced=already_balanced, oe_normalize=oe_normalize,
        save_path=save_path, chrom=chrom,
    )

    force = build_boltzmann_force(
        H, rc,
        alpha=alpha,
        k_scale=k_scale,
        force_group=force_group,
        return_controller=return_controller,
        noise_seed=noise_seed,
        r_min=r_min,
        kernel=kernel,
    )

    log.info("build_hic_force: done (group %d)", force_group)
    log.info("═" * 60)
    return force

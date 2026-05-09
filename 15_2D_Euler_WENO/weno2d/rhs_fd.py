"""
Finite-Difference WENO RHS for the 2-D compressible Euler equations.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

Implements the finite-difference WENO method (Shu 2009): the physical
flux is split into positive and negative components F = F⁺ + F⁻, and
each component is reconstructed independently with a WENO polynomial to
obtain the numerical flux at each cell interface.

Supported flux splittings
--------------------------
'LF'  — global Lax–Friedrichs: F⁺ = (F + a_max * q) / 2
'LLF' — local  Lax–Friedrichs (Rusanov): a_max computed cell-by-cell

Key performance improvements over the original MATLAB code
----------------------------------------------------------
Same vectorisation strategy as rhs_fv.py: the face-by-face loops and the
component-wise ``parfor`` loops are replaced by a single NumPy operation.
"""

import numpy as np

from .reconstructions import get_recon_x, get_recon_y, stencil_radius
from .flux_solvers import lf_flux_splitting, rusanov_flux_splitting
from .boundary_conditions import apply_bcs


def fd_weno_rhs(q, a_global, nx, ny, dx, dy, t,
                split_method, recon_method, test,
                gamma,
                mesh_wedge_position=None, preshock=None, postshock=None,
                shock_speed=None):
    """
    Compute the finite-difference WENO residual  L(q)  for the 2-D Euler equations.

    Parameters
    ----------
    q            : ndarray, shape (ny, nx, 4)
                   Conserved state including ghost cells; modified in place.
    a_global     : float
                   Global maximum wave speed (required for 'LF' splitting).
    nx, ny       : int — total cells in x and y (including ghost cells).
    dx, dy       : float — uniform cell widths.
    t            : float — current simulation time.
    split_method : str — 'LF' (global) or 'LLF' (local/Rusanov) flux splitting.
    recon_method : str — 'WENO5', 'WENO7', 'Poly5', or 'Poly7'.
    test         : str — 'Smooth', 'Riemann', or 'DMR'.
    gamma        : float — ratio of specific heats.
    mesh_wedge_position, preshock, postshock, shock_speed
                 : DMR-specific parameters.

    Returns
    -------
    res : ndarray, shape (ny, nx, 4) — zero in ghost cells.
    """
    R       = stencil_radius(recon_method)
    recon_x = get_recon_x(recon_method)
    recon_y = get_recon_y(recon_method)

    if split_method == 'LF':
        split_fn = lambda q_sub, normal: lf_flux_splitting(q_sub, normal, gamma, a_global)
    elif split_method == 'LLF':
        split_fn = lambda q_sub, normal: rusanov_flux_splitting(q_sub, normal, gamma)
    else:
        raise ValueError(f"Unknown flux splitting '{split_method}'. "
                         "Choose 'LF' or 'LLF'.")

    # 1. Apply boundary conditions
    apply_bcs(q, R, nx, ny, dx, dy, t, test,
              mesh_wedge_position=mesh_wedge_position,
              preshock=preshock, postshock=postshock,
              shock_speed=shock_speed)

    res = np.zeros_like(q)

    # ------------------------------------------------------------------ #
    # X-direction sweep                                                    #
    # ------------------------------------------------------------------ #
    # Flux splitting on interior y-rows, full x extent (includes ghosts)
    q_ic = q[R:ny-R, :, :]          # shape (nc_y, nx, 4)
    fp_x, fm_x = split_fn(q_ic, (1.0, 0.0))
    # fp_x, fm_x each shape (nc_y, nx, 4)

    # WENO reconstruction on each component independently
    flux_list = []
    for e in range(4):
        fl, fr = recon_x(fp_x[:, :, e], fm_x[:, :, e], nx)
        flux_list.append(fl + fr)        # combined numerical flux at each face
    # Each element shape (nc_y, nf_x); stack → (nc_y, nf_x, 4)
    flux_x = np.stack(flux_list, axis=-1)

    res[R:ny-R, R:nx-R, :] += (flux_x[:, :-1, :] - flux_x[:, 1:, :]) / dx

    # ------------------------------------------------------------------ #
    # Y-direction sweep                                                    #
    # ------------------------------------------------------------------ #
    q_jc = q[:, R:nx-R, :]          # shape (ny, nc_x, 4)
    fp_y, fm_y = split_fn(q_jc, (0.0, 1.0))

    flux_list = []
    for e in range(4):
        fl, fr = recon_y(fp_y[:, :, e], fm_y[:, :, e], ny)
        flux_list.append(fl + fr)
    flux_y = np.stack(flux_list, axis=-1)  # (nf_y, nc_x, 4)

    res[R:ny-R, R:nx-R, :] += (flux_y[:-1, :, :] - flux_y[1:, :, :]) / dy

    return res

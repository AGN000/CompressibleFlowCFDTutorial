"""
Finite-Volume WENO RHS for the 2-D compressible Euler equations.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

Implements the method-of-lines semi-discrete form

    dq/dt = -∂F/∂x - ∂G/∂y

using WENO reconstruction of interface states and an exact or approximate
Riemann solver at each interface.

Key performance improvements over the original MATLAB code
----------------------------------------------------------
1.  NumPy slice-based stencil extraction replaces reshape/flatten operations.
2.  Riemann solvers are evaluated on *all* interfaces simultaneously — the
    face-by-face loop ``for c in range(N)`` is eliminated.
3.  The residual accumulation is a single vectorised finite-difference
    operation: res[R:ny-R, R:nx-R] += (flux[:, :-1] - flux[:, 1:]) / dx.
4.  The ``parfor`` loop over the four conserved-variable components is
    replaced by a single call that processes all components in one pass.
"""

import numpy as np

from .reconstructions import get_recon_x, get_recon_y, stencil_radius
from .flux_solvers import get_flux_solver
from .boundary_conditions import apply_bcs


def fv_weno_rhs(q, a_global, nx, ny, dx, dy, t,
                flux_method, recon_method, test,
                gamma,
                mesh_wedge_position=None, preshock=None, postshock=None,
                shock_speed=None):
    """
    Compute the finite-volume WENO residual  L(q)  for the 2-D Euler equations.

    The residual is the right-hand side of  dq/dt = L(q), i.e. it represents
    -(∂F/∂x + ∂G/∂y) evaluated at cell centres.

    Parameters
    ----------
    q             : ndarray, shape (ny, nx, 4)
                    Conserved state including ghost cells; modified in place
                    to apply boundary conditions.
    a_global      : float
                    Global maximum wave speed (only used when flux_method='LF').
    nx, ny        : int — total cells in x and y (including ghost cells).
    dx, dy        : float — uniform cell widths.
    t             : float — current simulation time (for DMR top BC).
    flux_method   : str — 'LF', 'LLF', 'ROE', 'HLLE', or 'HLLC'.
    recon_method  : str — 'WENO5', 'WENO7', 'Poly5', or 'Poly7'.
    test          : str — 'Smooth', 'Riemann', or 'DMR'.
    gamma         : float — ratio of specific heats.
    mesh_wedge_position, preshock, postshock, shock_speed
                  : DMR-specific parameters (ignored for other tests).

    Returns
    -------
    res : ndarray, shape (ny, nx, 4) — the RHS residual (zero in ghost cells).
    """
    R     = stencil_radius(recon_method)
    recon_x = get_recon_x(recon_method)
    recon_y = get_recon_y(recon_method)
    flux_fn = get_flux_solver(flux_method)

    # 1. Apply boundary conditions (modifies q in place)
    apply_bcs(q, R, nx, ny, dx, dy, t, test,
              mesh_wedge_position=mesh_wedge_position,
              preshock=preshock, postshock=postshock,
              shock_speed=shock_speed)

    res = np.zeros_like(q)

    # ------------------------------------------------------------------ #
    # X-direction sweep                                                    #
    # ------------------------------------------------------------------ #
    # Interior y-rows only; full x-extent (ghost cells needed by stencil)
    q_ic = q[R:ny-R, :, :]          # shape (nc_y, nx, 4)  nc_y = ny-2R

    # Reconstruct all four components simultaneously; each output is
    # (nc_y, nf_x) where nf_x = nx - 2R + 1.
    qL_list, qR_list = [], []
    for e in range(4):
        w = q_ic[:, :, e]            # shape (nc_y, nx)
        qL_e, qR_e = recon_x(w, w, nx)
        qL_list.append(qL_e)
        qR_list.append(qR_e)

    # Stack → shape (nc_y, nf_x, 4)
    qL_x = np.stack(qL_list, axis=-1)
    qR_x = np.stack(qR_list, axis=-1)

    # Compute numerical flux at every x-interface in one vectorised call
    if flux_method == 'LF':
        flux_x = flux_fn(qL_x, qR_x, (1.0, 0.0), gamma, a_global)
    else:
        flux_x = flux_fn(qL_x, qR_x, (1.0, 0.0), gamma)
    # flux_x shape: (nc_y, nf_x, 4)

    # Divergence: res[j, i] += (flux[j,i] - flux[j,i+1]) / dx
    # Equivalent to the original double loop over cells and faces.
    res[R:ny-R, R:nx-R, :] += (flux_x[:, :-1, :] - flux_x[:, 1:, :]) / dx

    # ------------------------------------------------------------------ #
    # Y-direction sweep                                                    #
    # ------------------------------------------------------------------ #
    q_jc = q[:, R:nx-R, :]          # shape (ny, nc_x, 4)  nc_x = nx-2R

    qL_list, qR_list = [], []
    for e in range(4):
        w = q_jc[:, :, e]            # shape (ny, nc_x)
        qL_e, qR_e = recon_y(w, w, ny)
        qL_list.append(qL_e)
        qR_list.append(qR_e)

    # Stack → shape (nf_y, nc_x, 4)  nf_y = ny - 2R + 1
    qL_y = np.stack(qL_list, axis=-1)
    qR_y = np.stack(qR_list, axis=-1)

    if flux_method == 'LF':
        flux_y = flux_fn(qL_y, qR_y, (0.0, 1.0), gamma, a_global)
    else:
        flux_y = flux_fn(qL_y, qR_y, (0.0, 1.0), gamma)
    # flux_y shape: (nf_y, nc_x, 4)

    res[R:ny-R, R:nx-R, :] += (flux_y[:-1, :, :] - flux_y[1:, :, :]) / dy

    return res

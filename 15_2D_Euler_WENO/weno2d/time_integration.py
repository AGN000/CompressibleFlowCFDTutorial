"""
Time integration for the 2-D WENO Euler solver.

Provides the 3-stage, 3rd-order Strong-Stability-Preserving Runge–Kutta
method (SSP-RK3 / Shu–Osher) used in all driver scripts.

Reference
---------
Shu, C.-W. and Osher, S., "Efficient implementation of essentially
non-oscillatory shock-capturing schemes," J. Comput. Phys. 77, 439–471, 1988.
"""

import numpy as np


def ssprk3_step(q, dt, rhs_fn):
    """
    Advance the state *q* by one time step *dt* using SSP-RK3.

    The three-stage update is:

        q¹  = qⁿ + dt · L(qⁿ)
        q²  = ¾ qⁿ + ¼ (q¹ + dt · L(q¹))
        qⁿ⁺¹ = ⅓ qⁿ + ⅔ (q² + dt · L(q²))

    This is the convex-combination form that preserves the TVD/TVB property.

    Parameters
    ----------
    q      : ndarray, shape (ny, nx, 4)
             Conserved state at time level n (including ghost cells).
             A copy is made internally; the input is not modified.
    dt     : float — time step size.
    rhs_fn : callable — ``rhs_fn(q) → residual`` of same shape as q.
             The boundary conditions are applied inside rhs_fn.

    Returns
    -------
    q_new : ndarray, same shape as q — state at time level n+1.
    """
    q0 = q.copy()

    # Stage 1
    L1 = rhs_fn(q0.copy())
    q1 = q0 + dt*L1

    # Stage 2
    L2 = rhs_fn(q1)
    q2 = 0.75*q0 + 0.25*(q1 + dt*L2)

    # Stage 3
    L3 = rhs_fn(q2)
    return (q0 + 2.0*(q2 + dt*L3)) / 3.0


def cfl_dt(q, dx, dy, gamma, cfl, R):
    """
    Compute the maximum stable time step from the CFL condition.

    Parameters
    ----------
    q     : ndarray, shape (ny, nx, 4) — conserved state.
    dx,dy : float — cell widths.
    gamma : float — ratio of specific heats.
    cfl   : float — CFL number (typically 0.4 – 0.6 for SSP-RK3 + WENO5).
    R     : int   — ghost-cell radius (used to restrict to interior cells).

    Returns
    -------
    dt      : float — time step.
    a_global: float — global maximum wave speed (for LF flux splitting).
    """
    ny, nx = q.shape[:2]
    q_int = q[R:ny-R, R:nx-R, :]
    r = q_int[..., 0]
    u = q_int[..., 1] / r
    v = q_int[..., 2] / r
    E = q_int[..., 3]
    p = (gamma - 1.0)*(E - 0.5*r*(u**2 + v**2))
    c = np.sqrt(np.maximum(gamma*p/r, 0.0))

    spd = np.sqrt(u**2 + v**2) + c
    a_global = float(np.max(spd))
    dt = cfl * min(dx, dy) / a_global
    return dt, a_global

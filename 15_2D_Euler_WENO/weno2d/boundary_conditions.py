"""
Ghost-cell boundary conditions for the 2-D WENO Euler solver.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

All routines modify the conserved-variable array *q* (shape (ny, nx, 4))
in place.  Interior cells occupy indices [R:ny-R, R:nx-R, :].

Supported test cases
--------------------
'Smooth'  — fully periodic (smooth vortex, convergence studies).
'Riemann' — zero-order (Neumann) outflow on all four boundaries.
'DMR'     — Double Mach Reflection: Neumann left/right, reflective wall at
            bottom beyond wedge position, time-dependent moving shock at top.
"""

import numpy as np


def _distance_to_shock(x, y, t, shock_speed):
    """
    Signed distance from point (x, y) to the oblique shock front at time t.

    The shock travels at angle 30° to the x-axis (tan = 1/sqrt(3)) and
    starts at x = 1/6.  Positive = ahead (undisturbed), negative = behind.
    """
    slope = 1.0 / np.tan(np.pi / 6.0)   # = sqrt(3)
    x0    = 1.0 / 6.0
    return (slope*(x - x0 - shock_speed*t) - y) / np.sqrt(slope**2 + 1.0)


def apply_bcs(q, R, nx, ny, dx, dy, t, test,
              mesh_wedge_position=None, preshock=None, postshock=None,
              shock_speed=None):
    """
    Fill ghost-cell layers of *q* according to the boundary condition *test*.

    Parameters
    ----------
    q                  : ndarray (ny, nx, 4), modified in place.
    R                  : int — ghost-cell radius (= stencil radius).
    nx, ny             : int — total cells including ghost cells.
    dx, dy             : float — cell widths.
    t                  : float — current simulation time.
    test               : str — 'Smooth', 'Riemann', or 'DMR'.
    mesh_wedge_position: float — wedge x-position in cell units (DMR only).
    preshock           : ndarray (4,) — far-field state ahead of shock (DMR).
    postshock          : ndarray (4,) — state behind shock (DMR).
    shock_speed        : float — shock propagation speed in x (DMR).
    """
    if test == 'Smooth':
        # Periodic in x
        q[:, :R,      :] = q[:, nx-2*R:nx-R,   :]
        q[:, nx-R:,   :] = q[:, R:2*R,          :]
        # Periodic in y
        q[:R,      :, :] = q[ny-2*R:ny-R, :,    :]
        q[ny-R:,   :, :] = q[R:2*R,       :,    :]

    elif test == 'Riemann':
        # Neumann (zero-gradient) outflow on all sides
        q[:, :R,      :] = q[:, [R],       :]   # left  ghost ← first interior col
        q[:, nx-R:,   :] = q[:, [nx-R-1],  :]   # right ghost ← last  interior col
        q[:R,      :, :] = q[[R],       :, :]   # bottom ghost ← first interior row
        q[ny-R:,   :, :] = q[[ny-R-1],  :, :]   # top ghost   ← last  interior row

    elif test == 'DMR':
        # Left and right: Neumann
        q[:, :R,    :] = q[:, [R],      :]
        q[:, nx-R:, :] = q[:, [nx-R-1], :]

        # Bottom ghost rows — reflective wall beyond the wedge
        # Copy the first interior row into all bottom ghost layers
        q[:R, R:nx-R, :] = q[[R], R:nx-R, :]
        # Negate v-momentum (index 2) beyond wedge: reflective BC
        i_wedge = R + int(mesh_wedge_position)   # 0-indexed x threshold
        q[:R, i_wedge:nx-R, 2] = -q[[R], i_wedge:nx-R, 2]

        # Top ghost rows — full Dirichlet based on which side of the moving shock
        # each ghost cell lies on.  Pre-shock cells get preshock, post-shock
        # cells get postshock.  This is the correct BC for DMR (Woodward & Colella).
        j_top = np.arange(ny-R, ny)             # shape (R,)
        i_int = np.arange(R, nx-R)              # shape (nx-2R,)
        # Cell-centre coordinates (FV convention: offset by 0.5 cell from origin)
        x_pts = (i_int - R + 0.5) * dx          # x-coords of interior cell centres
        y_pts = (j_top[:, np.newaxis] - R + 0.5) * dy  # y-coords of ghost cell centres

        dist = _distance_to_shock(x_pts[np.newaxis, :], y_pts, t, shock_speed)

        for k in range(R):
            jj = ny - R + k            # actual row index in q
            row_dist = dist[k]         # shape (nx-2R,)
            post_mask = row_dist < 0.0  # behind shock → post-shock state
            q[jj, R:nx-R, :][post_mask]  = postshock
            q[jj, R:nx-R, :][~post_mask] = preshock

    else:
        raise ValueError(f"Unknown test case '{test}'. "
                         "Choose from 'Smooth', 'Riemann', 'DMR'.")

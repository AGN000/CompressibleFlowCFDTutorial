"""
WENO and polynomial interface reconstructions for the 2-D Euler solver.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

All functions operate on 2-D NumPy arrays and are fully vectorised — no
Python loops over individual cells.  The stencil notation follows Shu (2009):

    i-1/2        i+1/2
      |  I{i-1}  |  I{i}  |  I{i+1}  |
                 ↑
            face location

Left-biased  (qL) reconstruction uses cells shifted to the left.
Right-biased (qR) reconstruction uses cells shifted to the right.

Array shape conventions
-----------------------
X-direction routines  : input  (nr, N)   → output  (nr, nf)
Y-direction routines  : input  (N,  nc)  → output  (nf, nc)

where
    N   = total number of cells along the sweep direction (incl. ghost cells)
    nr  = number of rows    (for X sweep: interior rows only)
    nc  = number of columns (for Y sweep: interior columns only)
    nf  = number of interfaces = N - 2*R + 1
          (WENO5/Poly5: R=3 → nf = N-5;  WENO7/Poly7: R=4 → nf = N-7)

For the FD solver the same kernels accept the split-flux arrays fp and fm in
place of the conserved-variable array; the combined numerical flux is returned.
"""

import numpy as np

_EPS = 1e-6  # smoothness-indicator regularisation


# ---------------------------------------------------------------------------
# WENO5 — 5th-order weighted essentially non-oscillatory (R = 3)
# ---------------------------------------------------------------------------

def weno5_recon_x(v_left, v_right, N):
    """
    WENO5 reconstruction at all x-interfaces for a 2-D field.

    For FV use: pass the same state array for both v_left and v_right.
    For FD use: pass fp (positive split flux) for v_left and fm (negative
    split flux) for v_right; the return value is the combined numerical flux.

    Parameters
    ----------
    v_left  : ndarray, shape (nr, N)
        Field values used to build the left-biased (qL) reconstruction.
    v_right : ndarray, shape (nr, N)
        Field values used to build the right-biased (qR) reconstruction.
    N       : int
        Total number of cells in x (including ghost cells).

    Returns
    -------
    qL : ndarray, shape (nr, N-5)
        Left-biased interface values  q_{i+1/2}^{-}.
    qR : ndarray, shape (nr, N-5)
        Right-biased interface values q_{i+1/2}^{+}.
    """
    # Left-biased stencil  (uses cells i-2 … i+2)
    vmm = v_left[:, 0:N-5]
    vm  = v_left[:, 1:N-4]
    vo  = v_left[:, 2:N-3]
    vp  = v_left[:, 3:N-2]
    vpp = v_left[:, 4:N-1]

    B0n = 13/12*(vmm - 2*vm  + vo )**2 + 1/4*(vmm - 4*vm  + 3*vo)**2
    B1n = 13/12*(vm  - 2*vo  + vp )**2 + 1/4*(vm  - vp        )**2
    B2n = 13/12*(vo  - 2*vp  + vpp)**2 + 1/4*(3*vo - 4*vp + vpp)**2

    a0n = (1/10) / (_EPS + B0n)**2
    a1n = (6/10) / (_EPS + B1n)**2
    a2n = (3/10) / (_EPS + B2n)**2
    asum = a0n + a1n + a2n

    qL = (a0n*(2*vmm - 7*vm  + 11*vo) +
          a1n*( -vm  + 5*vo  +  2*vp) +
          a2n*(2*vo  + 5*vp  -    vpp)) / (6*asum)

    # Right-biased stencil (uses cells i-1 … i+3)
    umm = v_right[:, 1:N-4]
    um  = v_right[:, 2:N-3]
    uo  = v_right[:, 3:N-2]
    up  = v_right[:, 4:N-1]
    upp = v_right[:, 5:N  ]

    B0p = 13/12*(umm - 2*um  + uo )**2 + 1/4*(umm - 4*um  + 3*uo)**2
    B1p = 13/12*(um  - 2*uo  + up )**2 + 1/4*(um  - up        )**2
    B2p = 13/12*(uo  - 2*up  + upp)**2 + 1/4*(3*uo - 4*up + upp)**2

    a0p = (3/10) / (_EPS + B0p)**2
    a1p = (6/10) / (_EPS + B1p)**2
    a2p = (1/10) / (_EPS + B2p)**2
    asum = a0p + a1p + a2p

    qR = (a0p*(-umm + 5*um  +  2*uo) +
          a1p*(2*um  + 5*uo  -   up ) +
          a2p*(11*uo - 7*up  + 2*upp)) / (6*asum)

    return qL, qR


def weno5_recon_y(v_left, v_right, N):
    """
    WENO5 reconstruction at all y-interfaces for a 2-D field.

    Parameters
    ----------
    v_left  : ndarray, shape (N, nc)
    v_right : ndarray, shape (N, nc)
    N       : int — total cells in y including ghosts.

    Returns
    -------
    qL, qR : each ndarray of shape (N-5, nc).
    """
    # Left-biased stencil
    vmm = v_left[0:N-5, :]
    vm  = v_left[1:N-4, :]
    vo  = v_left[2:N-3, :]
    vp  = v_left[3:N-2, :]
    vpp = v_left[4:N-1, :]

    B0n = 13/12*(vmm - 2*vm  + vo )**2 + 1/4*(vmm - 4*vm  + 3*vo)**2
    B1n = 13/12*(vm  - 2*vo  + vp )**2 + 1/4*(vm  - vp        )**2
    B2n = 13/12*(vo  - 2*vp  + vpp)**2 + 1/4*(3*vo - 4*vp + vpp)**2

    a0n = (1/10) / (_EPS + B0n)**2
    a1n = (6/10) / (_EPS + B1n)**2
    a2n = (3/10) / (_EPS + B2n)**2
    asum = a0n + a1n + a2n

    qL = (a0n*(2*vmm - 7*vm  + 11*vo) +
          a1n*( -vm  + 5*vo  +  2*vp) +
          a2n*(2*vo  + 5*vp  -    vpp)) / (6*asum)

    # Right-biased stencil
    umm = v_right[1:N-4, :]
    um  = v_right[2:N-3, :]
    uo  = v_right[3:N-2, :]
    up  = v_right[4:N-1, :]
    upp = v_right[5:N,   :]

    B0p = 13/12*(umm - 2*um  + uo )**2 + 1/4*(umm - 4*um  + 3*uo)**2
    B1p = 13/12*(um  - 2*uo  + up )**2 + 1/4*(um  - up        )**2
    B2p = 13/12*(uo  - 2*up  + upp)**2 + 1/4*(3*uo - 4*up + upp)**2

    a0p = (3/10) / (_EPS + B0p)**2
    a1p = (6/10) / (_EPS + B1p)**2
    a2p = (1/10) / (_EPS + B2p)**2
    asum = a0p + a1p + a2p

    qR = (a0p*(-umm + 5*um  +  2*uo) +
          a1p*(2*um  + 5*uo  -   up ) +
          a2p*(11*uo - 7*up  + 2*upp)) / (6*asum)

    return qL, qR


# ---------------------------------------------------------------------------
# WENO7 — 7th-order weighted essentially non-oscillatory (R = 4)
# ---------------------------------------------------------------------------

def weno7_recon_x(v_left, v_right, N):
    """
    WENO7 reconstruction at all x-interfaces.

    Parameters
    ----------
    v_left, v_right : ndarray, shape (nr, N)
    N               : int

    Returns
    -------
    qL, qR : each ndarray of shape (nr, N-7).
    """
    # Left-biased stencil (cells i-3 … i+3)
    vmmm = v_left[:, 0:N-7]
    vmm  = v_left[:, 1:N-6]
    vm   = v_left[:, 2:N-5]
    vo   = v_left[:, 3:N-4]
    vp   = v_left[:, 4:N-3]
    vpp  = v_left[:, 5:N-2]
    vppp = v_left[:, 6:N-1]

    B0n = (vm*(134241*vm - 114894*vo) +
           vmmm*(56694*vm - 47214*vmm + 6649*vmmm - 22778*vo) +
           25729*vo**2 + vmm*(-210282*vm + 85641*vmm + 86214*vo))
    B1n = (vo*(41001*vo - 30414*vp) +
           vmm*(-19374*vm + 3169*vmm + 19014*vo - 5978*vp) +
           6649*vp**2 + vm*(33441*vm - 70602*vo + 23094*vp))
    B2n = (vp*(33441*vp - 19374*vpp) +
           vm*(6649*vm - 30414*vo + 23094*vp - 5978*vpp) +
           3169*vpp**2 + vo*(41001*vo - 70602*vp + 19014*vpp))
    B3n = (vpp*(85641*vpp - 47214*vppp) +
           vo*(25729*vo - 114894*vp + 86214*vpp - 22778*vppp) +
           6649*vppp**2 + vp*(134241*vp - 210282*vpp + 56694*vppp))

    g0, g1, g2, g3 = 1/35, 12/35, 18/35, 4/35
    a0n = g0 / (_EPS + B0n)**2
    a1n = g1 / (_EPS + B1n)**2
    a2n = g2 / (_EPS + B2n)**2
    a3n = g3 / (_EPS + B3n)**2
    asum = a0n + a1n + a2n + a3n

    qL = (a0n*(-3*vmmm + 13*vmm - 23*vm  + 25*vo            ) +
          a1n*(   vmm  -  5*vm  + 13*vo  +  3*vp             ) +
          a2n*(       -    vm   +  7*vo  +  7*vp  -    vpp   ) +
          a3n*(              3*vo + 13*vp -  5*vpp +   vppp  )) / (12*asum)

    # Right-biased stencil (cells i-2 … i+4)
    ummm = v_right[:, 1:N-6]
    umm  = v_right[:, 2:N-5]
    um   = v_right[:, 3:N-4]
    uo   = v_right[:, 4:N-3]
    up   = v_right[:, 5:N-2]
    upp  = v_right[:, 6:N-1]
    uppp = v_right[:, 7:N  ]

    B0p = (um*(134241*um - 114894*uo) +
           ummm*(56694*um - 47214*umm + 6649*ummm - 22778*uo) +
           25729*uo**2 + umm*(-210282*um + 85641*umm + 86214*uo))
    B1p = (uo*(41001*uo - 30414*up) +
           umm*(-19374*um + 3169*umm + 19014*uo - 5978*up) +
           6649*up**2 + um*(33441*um - 70602*uo + 23094*up))
    B2p = (up*(33441*up - 19374*upp) +
           um*(6649*um - 30414*uo + 23094*up - 5978*upp) +
           3169*upp**2 + uo*(41001*uo - 70602*up + 19014*upp))
    B3p = (upp*(85641*upp - 47214*uppp) +
           uo*(25729*uo - 114894*up + 86214*upp - 22778*uppp) +
           6649*uppp**2 + up*(134241*up - 210282*upp + 56694*uppp))

    g0, g1, g2, g3 = 4/35, 18/35, 12/35, 1/35
    a0p = g0 / (_EPS + B0p)**2
    a1p = g1 / (_EPS + B1p)**2
    a2p = g2 / (_EPS + B2p)**2
    a3p = g3 / (_EPS + B3p)**2
    asum = a0p + a1p + a2p + a3p

    qR = (a0p*(   ummm -  5*umm + 13*um  +  3*uo             ) +
          a1p*(        -   umm  +  7*um  +  7*uo  -    up    ) +
          a2p*(            3*um + 13*uo  -  5*up  +    upp   ) +
          a3p*(            25*uo - 23*up + 13*upp -  3*uppp  )) / (12*asum)

    return qL, qR


def weno7_recon_y(v_left, v_right, N):
    """
    WENO7 reconstruction at all y-interfaces.

    Parameters
    ----------
    v_left, v_right : ndarray, shape (N, nc)
    N               : int

    Returns
    -------
    qL, qR : each ndarray of shape (N-7, nc).
    """
    vmmm = v_left[0:N-7, :]
    vmm  = v_left[1:N-6, :]
    vm   = v_left[2:N-5, :]
    vo   = v_left[3:N-4, :]
    vp   = v_left[4:N-3, :]
    vpp  = v_left[5:N-2, :]
    vppp = v_left[6:N-1, :]

    B0n = (vm*(134241*vm - 114894*vo) +
           vmmm*(56694*vm - 47214*vmm + 6649*vmmm - 22778*vo) +
           25729*vo**2 + vmm*(-210282*vm + 85641*vmm + 86214*vo))
    B1n = (vo*(41001*vo - 30414*vp) +
           vmm*(-19374*vm + 3169*vmm + 19014*vo - 5978*vp) +
           6649*vp**2 + vm*(33441*vm - 70602*vo + 23094*vp))
    B2n = (vp*(33441*vp - 19374*vpp) +
           vm*(6649*vm - 30414*vo + 23094*vp - 5978*vpp) +
           3169*vpp**2 + vo*(41001*vo - 70602*vp + 19014*vpp))
    B3n = (vpp*(85641*vpp - 47214*vppp) +
           vo*(25729*vo - 114894*vp + 86214*vpp - 22778*vppp) +
           6649*vppp**2 + vp*(134241*vp - 210282*vpp + 56694*vppp))

    g0, g1, g2, g3 = 1/35, 12/35, 18/35, 4/35
    a0n = g0 / (_EPS + B0n)**2
    a1n = g1 / (_EPS + B1n)**2
    a2n = g2 / (_EPS + B2n)**2
    a3n = g3 / (_EPS + B3n)**2
    asum = a0n + a1n + a2n + a3n

    qL = (a0n*(-3*vmmm + 13*vmm - 23*vm  + 25*vo             ) +
          a1n*(   vmm  -  5*vm  + 13*vo  +  3*vp              ) +
          a2n*(       -    vm   +  7*vo  +  7*vp  -    vpp    ) +
          a3n*(              3*vo + 13*vp -  5*vpp +   vppp   )) / (12*asum)

    ummm = v_right[1:N-6, :]
    umm  = v_right[2:N-5, :]
    um   = v_right[3:N-4, :]
    uo   = v_right[4:N-3, :]
    up   = v_right[5:N-2, :]
    upp  = v_right[6:N-1, :]
    uppp = v_right[7:N,   :]

    B0p = (um*(134241*um - 114894*uo) +
           ummm*(56694*um - 47214*umm + 6649*ummm - 22778*uo) +
           25729*uo**2 + umm*(-210282*um + 85641*umm + 86214*uo))
    B1p = (uo*(41001*uo - 30414*up) +
           umm*(-19374*um + 3169*umm + 19014*uo - 5978*up) +
           6649*up**2 + um*(33441*um - 70602*uo + 23094*up))
    B2p = (up*(33441*up - 19374*upp) +
           um*(6649*um - 30414*uo + 23094*up - 5978*upp) +
           3169*upp**2 + uo*(41001*uo - 70602*up + 19014*upp))
    B3p = (upp*(85641*upp - 47214*uppp) +
           uo*(25729*uo - 114894*up + 86214*upp - 22778*uppp) +
           6649*uppp**2 + up*(134241*up - 210282*upp + 56694*uppp))

    g0, g1, g2, g3 = 4/35, 18/35, 12/35, 1/35
    a0p = g0 / (_EPS + B0p)**2
    a1p = g1 / (_EPS + B1p)**2
    a2p = g2 / (_EPS + B2p)**2
    a3p = g3 / (_EPS + B3p)**2
    asum = a0p + a1p + a2p + a3p

    qR = (a0p*(   ummm -  5*umm + 13*um  +  3*uo              ) +
          a1p*(        -   umm  +  7*um  +  7*uo  -    up     ) +
          a2p*(            3*um + 13*uo  -  5*up  +    upp    ) +
          a3p*(            25*uo - 23*up + 13*upp -  3*uppp   )) / (12*asum)

    return qL, qR


# ---------------------------------------------------------------------------
# Polynomial reconstructions (no nonlinear weights — linear, fixed stencil)
# ---------------------------------------------------------------------------

def poly5_recon_x(v_left, v_right, N):
    """5th-order polynomial reconstruction in x (R=3 stencil)."""
    vmm = v_left[:, 0:N-5]; vm = v_left[:, 1:N-4]
    vo  = v_left[:, 2:N-3]; vp = v_left[:, 3:N-2]; vpp = v_left[:, 4:N-1]
    qL = (2*vmm - 13*vm + 47*vo + 27*vp - 3*vpp) / 60

    umm = v_right[:, 1:N-4]; um = v_right[:, 2:N-3]
    uo  = v_right[:, 3:N-2]; up = v_right[:, 4:N-1]; upp = v_right[:, 5:N]
    qR = (-3*umm + 27*um + 47*uo - 13*up + 2*upp) / 60
    return qL, qR


def poly5_recon_y(v_left, v_right, N):
    """5th-order polynomial reconstruction in y (R=3 stencil)."""
    vmm = v_left[0:N-5, :]; vm = v_left[1:N-4, :]
    vo  = v_left[2:N-3, :]; vp = v_left[3:N-2, :]; vpp = v_left[4:N-1, :]
    qL = (2*vmm - 13*vm + 47*vo + 27*vp - 3*vpp) / 60

    umm = v_right[1:N-4, :]; um = v_right[2:N-3, :]
    uo  = v_right[3:N-2, :]; up = v_right[4:N-1, :]; upp = v_right[5:N, :]
    qR = (-3*umm + 27*um + 47*uo - 13*up + 2*upp) / 60
    return qL, qR


def poly7_recon_x(v_left, v_right, N):
    """7th-order polynomial reconstruction in x (R=4 stencil)."""
    vmmm = v_left[:, 0:N-7]; vmm = v_left[:, 1:N-6]; vm = v_left[:, 2:N-5]
    vo   = v_left[:, 3:N-4]; vp  = v_left[:, 4:N-3]; vpp = v_left[:, 5:N-2]
    vppp = v_left[:, 6:N-1]
    qL = (-3*vmmm + 25*vmm - 101*vm + 319*vo + 214*vp - 38*vpp + 4*vppp) / 420

    ummm = v_right[:, 1:N-6]; umm = v_right[:, 2:N-5]; um = v_right[:, 3:N-4]
    uo   = v_right[:, 4:N-3]; up  = v_right[:, 5:N-2]; upp = v_right[:, 6:N-1]
    uppp = v_right[:, 7:N]
    qR = (4*ummm - 38*umm + 214*um + 319*uo - 101*up + 25*upp - 3*uppp) / 420
    return qL, qR


def poly7_recon_y(v_left, v_right, N):
    """7th-order polynomial reconstruction in y (R=4 stencil)."""
    vmmm = v_left[0:N-7, :]; vmm = v_left[1:N-6, :]; vm = v_left[2:N-5, :]
    vo   = v_left[3:N-4, :]; vp  = v_left[4:N-3, :]; vpp = v_left[5:N-2, :]
    vppp = v_left[6:N-1, :]
    qL = (-3*vmmm + 25*vmm - 101*vm + 319*vo + 214*vp - 38*vpp + 4*vppp) / 420

    ummm = v_right[1:N-6, :]; umm = v_right[2:N-5, :]; um = v_right[3:N-4, :]
    uo   = v_right[4:N-3, :]; up  = v_right[5:N-2, :]; upp = v_right[6:N-1, :]
    uppp = v_right[7:N, :]
    qR = (4*ummm - 38*umm + 214*um + 319*uo - 101*up + 25*upp - 3*uppp) / 420
    return qL, qR


# ---------------------------------------------------------------------------
# Dispatch helpers
# ---------------------------------------------------------------------------

_RECON_X = {
    'WENO5': weno5_recon_x, 'WENO7': weno7_recon_x,
    'Poly5': poly5_recon_x, 'Poly7': poly7_recon_x,
}
_RECON_Y = {
    'WENO5': weno5_recon_y, 'WENO7': weno7_recon_y,
    'Poly5': poly5_recon_y, 'Poly7': poly7_recon_y,
}
_STENCIL_R = {'WENO5': 3, 'WENO7': 4, 'Poly5': 3, 'Poly7': 4}


def get_recon_x(name):
    """Return the X-direction reconstruction function for *name*."""
    try:
        return _RECON_X[name]
    except KeyError:
        raise ValueError(f"Unknown reconstruction '{name}'. "
                         f"Choose from {list(_RECON_X)}")


def get_recon_y(name):
    """Return the Y-direction reconstruction function for *name*."""
    try:
        return _RECON_Y[name]
    except KeyError:
        raise ValueError(f"Unknown reconstruction '{name}'. "
                         f"Choose from {list(_RECON_Y)}")


def stencil_radius(name):
    """Return the ghost-cell radius R for reconstruction *name*."""
    return _STENCIL_R[name]

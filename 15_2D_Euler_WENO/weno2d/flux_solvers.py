"""
Vectorised numerical flux functions for the 2-D compressible Euler equations.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

All routines accept left and right interface states as arrays and compute
the numerical flux at every interface simultaneously, eliminating the
face-by-face Python loop present in the original MATLAB code.

Conservative state vector
--------------------------
q = [rho, rho*u, rho*v, E]ᵀ

Physical fluxes in normal direction n = (nx, ny)
-------------------------------------------------
F·n = [rho*vn,  rho*vn*u + p*nx,  rho*vn*v + p*ny,  vn*(E+p)]
where  vn = u*nx + v*ny

Array conventions
-----------------
qL, qR : shape (..., 4) — last axis is the conserved-variable index.
normal  : length-2 tuple (nx, ny).
gamma   : ratio of specific heats (scalar).

Returns flux of the same shape as qL.
"""

import numpy as np


def _primitives(q, gamma):
    """Extract primitive variables and total enthalpy from conserved state q."""
    r  = q[..., 0]
    u  = q[..., 1] / r
    v  = q[..., 2] / r
    E  = q[..., 3]
    p  = (gamma - 1.0) * (E - 0.5*r*(u**2 + v**2))
    H  = (E + p) / r
    return r, u, v, p, H


def _euler_flux_n(r, u, v, p, H, nx, ny):
    """Physical Euler flux in normal direction; returns shape (..., 4)."""
    vn = u*nx + v*ny
    return np.stack([r*vn,
                     r*vn*u + p*nx,
                     r*vn*v + p*ny,
                     r*vn*H], axis=-1)


# ---------------------------------------------------------------------------
# Global Lax–Friedrichs (requires a global wave-speed bound passed in)
# ---------------------------------------------------------------------------

def lf_flux(qL, qR, normal, gamma, a_global):
    """
    Global Lax–Friedrichs numerical flux.

    Parameters
    ----------
    qL, qR   : ndarray, shape (..., 4)
    normal   : (nx, ny)
    gamma    : float
    a_global : float — global maximum wave speed for the current time step.

    Returns
    -------
    flux : ndarray, same shape as qL.
    """
    nx, ny = normal
    rL, uL, vL, pL, HL = _primitives(qL, gamma)
    rR, uR, vR, pR, HR = _primitives(qR, gamma)
    FL = _euler_flux_n(rL, uL, vL, pL, HL, nx, ny)
    FR = _euler_flux_n(rR, uR, vR, pR, HR, nx, ny)
    return 0.5*(FL + FR + a_global*(qL - qR))


# ---------------------------------------------------------------------------
# Local Lax–Friedrichs / Rusanov
# ---------------------------------------------------------------------------

def rusanov_flux(qL, qR, normal, gamma):
    """
    Rusanov (local Lax–Friedrichs) numerical flux.

    Uses the Roe-averaged wave speed as the local dissipation coefficient.

    Parameters
    ----------
    qL, qR : ndarray, shape (..., 4)
    normal : (nx, ny)
    gamma  : float

    Returns
    -------
    flux : ndarray, same shape as qL.
    """
    nx, ny = normal
    rL, uL, vL, pL, HL = _primitives(qL, gamma)
    rR, uR, vR, pR, HR = _primitives(qR, gamma)

    # Roe-averaged speed of sound for local wave-speed estimate
    RT  = np.sqrt(np.maximum(rR / rL, 0.0))
    u_r = (uL + RT*uR) / (1.0 + RT)
    v_r = (vL + RT*vR) / (1.0 + RT)
    H_r = (HL + RT*HR) / (1.0 + RT)
    a_r = np.sqrt(np.maximum((gamma - 1.0)*(H_r - 0.5*(u_r**2 + v_r**2)), 0.0))
    smax = np.abs(u_r*nx + v_r*ny) + a_r

    FL = _euler_flux_n(rL, uL, vL, pL, HL, nx, ny)
    FR = _euler_flux_n(rR, uR, vR, pR, HR, nx, ny)
    return 0.5*(FL + FR + smax[..., np.newaxis]*(qL - qR))


# ---------------------------------------------------------------------------
# Roe flux with Harten entropy fix
# ---------------------------------------------------------------------------

def roe_flux(qL, qR, normal, gamma):
    """
    Roe approximate Riemann solver with Harten's entropy fix.

    References
    ----------
    Roe, P. L., "Approximate Riemann solvers, parameter vectors, and difference
    schemes," J. Comput. Phys. 43, 357–372, 1981.
    Harten, A., "On a class of high resolution total-variation-stable finite
    difference schemes," SIAM J. Numer. Anal. 21, 1–23, 1984.

    Parameters
    ----------
    qL, qR : ndarray, shape (..., 4)
    normal : (nx, ny)
    gamma  : float

    Returns
    -------
    flux : ndarray, same shape as qL.
    """
    nx, ny = normal
    tx, ty = -ny, nx  # tangent direction

    rL, uL, vL, pL, HL = _primitives(qL, gamma)
    rR, uR, vR, pR, HR = _primitives(qR, gamma)

    vnL = uL*nx + vL*ny;  vtL = uL*tx + vL*ty
    vnR = uR*nx + vR*ny;  vtR = uR*tx + vR*ty

    RT  = np.sqrt(np.maximum(rR / rL, 0.0))
    r_r = RT * rL
    u_r = (uL + RT*uR) / (1.0 + RT)
    v_r = (vL + RT*vR) / (1.0 + RT)
    H_r = (HL + RT*HR) / (1.0 + RT)
    a_r = np.sqrt(np.maximum((gamma - 1.0)*(H_r - 0.5*(u_r**2 + v_r**2)), 0.0))
    vn_r = u_r*nx + v_r*ny
    vt_r = u_r*tx + v_r*ty

    # Jumps
    dr  = rR - rL
    dp  = pR - pL
    dvn = vnR - vnL
    dvt = vtR - vtL

    # Characteristic wave strengths
    dV0 = (dp - r_r*a_r*dvn) / (2.0*a_r**2)
    dV1 = r_r*dvt / a_r
    dV2 = dr - dp / a_r**2
    dV3 = (dp + r_r*a_r*dvn) / (2.0*a_r**2)

    # Wave speeds and Harten's entropy fix for nonlinear fields (λ1, λ4)
    ws0 = np.abs(vn_r - a_r)
    ws1 = np.abs(vn_r)
    ws3 = np.abs(vn_r + a_r)
    dws = 0.2  # entropy-fix threshold
    ws0 = np.where(ws0 < dws, 0.5*(ws0**2/dws + dws), ws0)
    ws3 = np.where(ws3 < dws, 0.5*(ws3**2/dws + dws), ws3)

    # Right eigenvectors (columns of R matrix)  — note: row 3 uses v for the
    # contact wave (original MATLAB code had a typo 'u' there).
    R0 = np.stack([np.ones_like(u_r), u_r - a_r*nx, v_r - a_r*ny, H_r - vn_r*a_r], axis=-1)
    R1 = np.stack([np.zeros_like(u_r), a_r*tx*np.ones_like(u_r), a_r*ty*np.ones_like(v_r), vt_r*a_r], axis=-1)
    R2 = np.stack([np.ones_like(u_r), u_r, v_r, 0.5*(u_r**2 + v_r**2)], axis=-1)
    R3 = np.stack([np.ones_like(u_r), u_r + a_r*nx, v_r + a_r*ny, H_r + vn_r*a_r], axis=-1)

    diss = (ws0[..., np.newaxis]*dV0[..., np.newaxis]*R0 +
            ws1[..., np.newaxis]*dV1[..., np.newaxis]*R1 +
            ws1[..., np.newaxis]*dV2[..., np.newaxis]*R2 +
            ws3[..., np.newaxis]*dV3[..., np.newaxis]*R3)

    FL = _euler_flux_n(rL, uL, vL, pL, HL, nx, ny)
    FR = _euler_flux_n(rR, uR, vR, pR, HR, nx, ny)
    return 0.5*(FL + FR - diss)


# ---------------------------------------------------------------------------
# HLLE (Einfeldt)
# ---------------------------------------------------------------------------

def hlle_flux(qL, qR, normal, gamma):
    """
    HLLE (Einfeldt) flux.

    Reference
    ---------
    Einfeldt, B., "On Godunov-type methods for gas dynamics,"
    SIAM J. Numer. Anal. 25(2), 294–318, 1988.

    Parameters
    ----------
    qL, qR : ndarray, shape (..., 4)
    normal : (nx, ny)
    gamma  : float

    Returns
    -------
    flux : ndarray, same shape as qL.
    """
    nx, ny = normal
    rL, uL, vL, pL, HL = _primitives(qL, gamma)
    rR, uR, vR, pR, HR = _primitives(qR, gamma)
    vnL = uL*nx + vL*ny
    vnR = uR*nx + vR*ny
    aL  = np.sqrt(np.maximum(gamma*pL/rL, 0.0))
    aR  = np.sqrt(np.maximum(gamma*pR/rR, 0.0))

    RT  = np.sqrt(np.maximum(rR / rL, 0.0))
    u_r = (uL + RT*uR) / (1.0 + RT)
    v_r = (vL + RT*vR) / (1.0 + RT)
    H_r = (HL + RT*HR) / (1.0 + RT)
    a_r = np.sqrt(np.maximum((gamma - 1.0)*(H_r - 0.5*(u_r**2 + v_r**2)), 0.0))
    vn_r = u_r*nx + v_r*ny

    SL = np.minimum(vnL - aL, vn_r - a_r)
    SR = np.maximum(vnR + aR, vn_r + a_r)
    SL = np.minimum(SL, 0.0)
    SR = np.maximum(SR, 0.0)

    FL = _euler_flux_n(rL, uL, vL, pL, HL, nx, ny)
    FR = _euler_flux_n(rR, uR, vR, pR, HR, nx, ny)
    denom = SR - SL
    return (SR[..., np.newaxis]*FL - SL[..., np.newaxis]*FR +
            (SL*SR)[..., np.newaxis]*(qR - qL)) / denom[..., np.newaxis]


# ---------------------------------------------------------------------------
# HLLC
# ---------------------------------------------------------------------------

def hllc_flux(qL, qR, normal, gamma):
    """
    HLLC flux (Toro 1994).

    Reference
    ---------
    Toro, E. F., Spruce, M., Speares, W., "Restoration of the contact surface
    in the HLL-Riemann solver," Shock Waves 4, 25–34, 1994.

    Parameters
    ----------
    qL, qR : ndarray, shape (..., 4)
    normal : (nx, ny)
    gamma  : float

    Returns
    -------
    flux : ndarray, same shape as qL.
    """
    nx, ny = normal
    rL, uL, vL, pL, HL = _primitives(qL, gamma)
    rR, uR, vR, pR, HR = _primitives(qR, gamma)
    vnL = uL*nx + vL*ny
    vnR = uR*nx + vR*ny

    FL = _euler_flux_n(rL, uL, vL, pL, HL, nx, ny)
    FR = _euler_flux_n(rR, uR, vR, pR, HR, nx, ny)

    # Clamp pressure/density for wave-speed estimates only — prevents NaN when
    # WENO5 produces slight oscillations that overshoot into negative pressure.
    _eps = 1e-14
    pLs = np.maximum(pL, _eps);  rLs = np.maximum(rL, _eps)
    pRs = np.maximum(pR, _eps);  rRs = np.maximum(rR, _eps)
    aL  = np.sqrt(gamma*pLs/rLs)
    aR  = np.sqrt(gamma*pRs/rRs)

    # Pressure estimate via PVRS (Toro)
    PPV  = np.maximum(_eps, 0.5*(pLs+pRs) + 0.5*(vnL-vnR)*0.25*(rLs+rRs)*(aL+aR))
    pmin = np.minimum(pLs, pRs)
    pmax = np.maximum(pLs, pRs)
    Qmax = pmax / pmin
    Quser = 2.0

    # Two-Rarefaction estimate
    PQ  = (pLs/pRs)**((gamma - 1.0)/(2.0*gamma))
    uTR = (PQ*vnL/aL + vnR/aR + 2.0/(gamma-1.0)*(PQ - 1.0)) / (PQ/aL + 1.0/aR)
    PTL = 1.0 + (gamma-1.0)/2.0*(vnL - uTR)/aL
    PTR = 1.0 + (gamma-1.0)/2.0*(uTR - vnR)/aR
    pTR = 0.5*(pLs*np.maximum(PTL, 0.0)**(2.0*gamma/(gamma-1.0)) +
               pRs*np.maximum(PTR, 0.0)**(2.0*gamma/(gamma-1.0)))

    # Two-Shock estimate
    GEL = np.sqrt((2.0/(gamma+1.0)/rLs) / ((gamma-1.0)/(gamma+1.0)*pLs + PPV))
    GER = np.sqrt((2.0/(gamma+1.0)/rRs) / ((gamma-1.0)/(gamma+1.0)*pRs + PPV))
    pTS = (GEL*pLs + GER*pRs - (vnR - vnL)) / (GEL + GER)

    use_PVRS = (Qmax <= Quser) & (pmin <= PPV) & (PPV <= pmax)
    use_TR   = (~use_PVRS) & (PPV < pmin)
    pM = np.where(use_PVRS, PPV, np.where(use_TR, pTR, pTS))

    # Wave speeds
    zL = np.where(pM > pLs, np.sqrt(1.0 + (gamma+1.0)/(2.0*gamma)*(pM/pLs - 1.0)), 1.0)
    zR = np.where(pM > pRs, np.sqrt(1.0 + (gamma+1.0)/(2.0*gamma)*(pM/pRs - 1.0)), 1.0)
    SL = vnL - aL*zL
    SR = vnR + aR*zR
    _denom_SM = rR*(SR - vnR) - rL*(SL - vnL)
    SM = (pL - pR + rR*vnR*(SR - vnR) - rL*vnL*(SL - vnL)) / \
         np.where(np.abs(_denom_SM) > 1e-14, _denom_SM, 1e-14)

    # Star states (use actual pL, pR — not clamped — for conservation)
    absny = np.abs(ny); absnx = np.abs(nx)
    _dSL_L = np.where(np.abs(SL - vnL) > 1e-14, SL - vnL, 1e-14)
    _dSR_R = np.where(np.abs(SR - vnR) > 1e-14, SR - vnR, 1e-14)
    facL = rL*(SL - vnL) / (SL - SM)
    qsL = np.stack([facL,
                    facL*(SM*nx + uL*absny),
                    facL*(SM*ny + vL*absnx),
                    facL*(qL[..., 3]/rL + (SM - vnL)*(SM + pL/(rL*_dSL_L)))],
                   axis=-1)
    facR = rR*(SR - vnR) / (SR - SM)
    qsR = np.stack([facR,
                    facR*(SM*nx + uR*absny),
                    facR*(SM*ny + vR*absnx),
                    facR*(qR[..., 3]/rR + (SM - vnR)*(SM + pR/(rR*_dSR_R)))],
                   axis=-1)

    # Select region
    FL_star = FL + SL[..., np.newaxis]*(qsL - qL)
    FR_star = FR + SR[..., np.newaxis]*(qsR - qR)

    flux = np.where((SL >= 0)[..., np.newaxis], FL,
           np.where((SM >= 0)[..., np.newaxis], FL_star,
           np.where((SR >= 0)[..., np.newaxis], FR_star, FR)))
    return flux


# ---------------------------------------------------------------------------
# Lax–Friedrichs flux splitting (for FD-WENO)
# ---------------------------------------------------------------------------

def lf_flux_splitting(q, normal, gamma, a_global):
    """
    Global Lax–Friedrichs flux splitting F = F⁺ + F⁻.

    Used by the finite-difference WENO scheme.  Returns fp and fm such that
    fp = (F + a*q)/2  and  fm = (F - a*q)/2,
    where F is the physical flux in direction *normal*.

    Parameters
    ----------
    q        : ndarray, shape (..., 4)
    normal   : (nx_normal, ny_normal)
    gamma    : float
    a_global : float

    Returns
    -------
    fp, fm : each ndarray of same shape as q.
    """
    nx, ny = normal
    r, u, v, p, H = _primitives(q, gamma)
    F = _euler_flux_n(r, u, v, p, H, nx, ny)
    fp = 0.5*(F + a_global*q)
    fm = 0.5*(F - a_global*q)
    return fp, fm


def rusanov_flux_splitting(q, normal, gamma):
    """
    Local Lax–Friedrichs (Rusanov) flux splitting F = F⁺ + F⁻.

    The local wave speed a(x) = |vn| + c is computed cell-by-cell.

    Parameters
    ----------
    q      : ndarray, shape (..., 4)
    normal : (nx_normal, ny_normal)
    gamma  : float

    Returns
    -------
    fp, fm : each ndarray of same shape as q.
    """
    nx, ny = normal
    r, u, v, p, H = _primitives(q, gamma)
    vn = u*nx + v*ny
    c  = np.sqrt(np.maximum(gamma*p/r, 0.0))
    a  = np.abs(vn) + c
    F  = _euler_flux_n(r, u, v, p, H, nx, ny)
    fp = 0.5*(F + a[..., np.newaxis]*q)
    fm = 0.5*(F - a[..., np.newaxis]*q)
    return fp, fm


# ---------------------------------------------------------------------------
# Dispatch helper
# ---------------------------------------------------------------------------

def get_flux_solver(name):
    """
    Return a flux-function callable for the given *name*.

    The returned function has signature  f(qL, qR, normal, gamma, **kw).
    For 'LF' an additional keyword *a_global* must be provided.
    """
    _map = {
        'LF':   lf_flux,
        'LLF':  rusanov_flux,
        'ROE':  roe_flux,
        'HLLE': hlle_flux,
        'HLLC': hllc_flux,
    }
    try:
        return _map[name]
    except KeyError:
        raise ValueError(f"Unknown flux method '{name}'. Choose from {list(_map)}")

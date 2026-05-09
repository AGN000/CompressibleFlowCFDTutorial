"""
weno2d — 2-D WENO/FV/FD solver for the compressible Euler equations.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO) with the
following changes:
  - Full NumPy vectorisation (no element-wise Python loops over cells/faces)
  - Unified WENO5/7 and Poly5/7 reconstruction kernels for X and Y directions
  - Vectorised Riemann solvers (LF, Rusanov, Roe, HLLE, HLLC)
  - SSP-RK3 time integrator decoupled from the RHS routines
  - Global variables replaced by explicit parameter passing

References
----------
[1] Shu, C.-W., "High order weighted essentially non-oscillatory schemes for
    convection dominated problems," SIAM Review 51(1), 82–126, 2009.
[2] Toro, E. F., "Riemann Solvers and Numerical Methods for Fluid Dynamics,"
    Springer, 2nd ed., 1999.
[3] Balsara, D. S., "A two-dimensional HLLC Riemann solver for conservation
    laws," J. Comput. Phys. 231, 7476–7503, 2012.
[4] Einfeldt, B., "On Godunov-type methods for gas dynamics," SIAM J. Numer.
    Anal. 25(2), 294–318, 1988.
"""

from .reconstructions import (
    weno5_recon_x, weno5_recon_y,
    weno7_recon_x, weno7_recon_y,
    poly5_recon_x, poly5_recon_y,
    poly7_recon_x, poly7_recon_y,
)
from .flux_solvers import lf_flux, rusanov_flux, roe_flux, hlle_flux, hllc_flux
from .rhs_fv import fv_weno_rhs
from .rhs_fd import fd_weno_rhs
from .initial_conditions import ic_double_mach_reflection, ic_riemann2d
from .time_integration import ssprk3_step

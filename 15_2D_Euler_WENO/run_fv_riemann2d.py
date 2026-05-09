"""
FV-WENO solver for 2-D Riemann problems.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

Solves the 2-D compressible Euler equations on the unit square using a
Finite-Volume WENO scheme with a user-selected Riemann solver.  Any of the
19 standard Kurganov–Tadmor configurations can be selected via the IC
parameter below.

Usage
-----
    python run_fv_riemann2d.py

Edit the "Parameters" block to change the resolution, scheme, or test case.

References
----------
[1] Kurganov, A. and Tadmor, E., "Solution of two-dimensional Riemann
    problems for gas dynamics without Riemann problem solvers,"
    Numer. Methods PDE 18(5), 584–608, 2002.
[2] Toro, E. F., "Riemann Solvers and Numerical Methods for Fluid Dynamics,"
    Springer, 2nd ed., 1999.
"""

import time
import numpy as np
import matplotlib
matplotlib.use('Agg')   # headless rendering; change to 'TkAgg' for interactive
import matplotlib.pyplot as plt

from weno2d import fv_weno_rhs, ic_riemann2d
from weno2d.reconstructions import stencil_radius
from weno2d.time_integration import ssprk3_step, cfl_dt

# ============================================================
# Parameters
# ============================================================
CFL        = 0.50          # CFL number
T_END      = 0.30          # final simulation time
NX         = 200           # physical cells in x
NY         = 200           # physical cells in y
N_DOF      = 5             # degrees of freedom (ideal diatomic = 5, monatomic = 3)
IC         = 6             # Riemann configuration (1 – 19, 'Sod_x', 'Sod_y')
FLUX_MTH   = 'HLLC'        # 'LF', 'LLF', 'ROE', 'HLLE', 'HLLC'
RECON_MTH  = 'WENO5'       # 'WENO5', 'WENO7', 'Poly5', 'Poly7'
SAVE_FIG   = True          # save final figure to disk
FIG_PATH   = 'riemann2d_result.png'

# ============================================================
# Derived physical parameters
# ============================================================
gamma = (N_DOF + 2.0) / N_DOF   # ratio of specific heats

# ============================================================
# Grid setup
# ============================================================
Lx = 1.0;  dx = Lx / NX
Ly = 1.0;  dy = Ly / NY
xc = np.linspace(dx/2.0, Lx - dx/2.0, NX)
yc = np.linspace(dy/2.0, Ly - dy/2.0, NY)
x, y = np.meshgrid(xc, yc)

# ============================================================
# Initial condition
# ============================================================
r0, u0, v0, p0 = ic_riemann2d(x, y, IC, gamma=gamma)
E0 = p0/(gamma - 1.0) + 0.5*r0*(u0**2 + v0**2)
c0 = np.sqrt(gamma*p0/r0)
Q0 = np.stack([r0, r0*u0, r0*v0, E0], axis=-1)   # shape (NY, NX, 4)

# ============================================================
# Ghost-cell padding
# ============================================================
R  = stencil_radius(RECON_MTH)
nx = NX + 2*R;  ny = NY + 2*R   # total cells including ghosts
q  = np.zeros((ny, nx, 4))
q[R:ny-R, R:nx-R, :] = Q0       # embed physical domain

# ============================================================
# Initial time step
# ============================================================
dt, a = cfl_dt(q, dx, dy, gamma, CFL, R)

# ============================================================
# Bind RHS with fixed parameters
# ============================================================
def rhs(q_in):
    return fv_weno_rhs(
        q_in, a, nx, ny, dx, dy, t,
        flux_method=FLUX_MTH,
        recon_method=RECON_MTH,
        test='Riemann',
        gamma=gamma,
    )

# ============================================================
# Time loop
# ============================================================
t     = 0.0
it    = 0
t_cpu = time.perf_counter()

print(f"\nFV-WENO Euler  |  Config {IC}  |  {RECON_MTH}-{FLUX_MTH}  |  "
      f"{NX}×{NY} cells  |  t_end={T_END}")
print(f"{'iter':>6s}  {'t':>10s}  {'dt':>10s}")

while t < T_END:
    dt = min(dt, T_END - t)
    q  = ssprk3_step(q, dt, rhs)
    t += dt

    # Recompute wave speed and time step after update
    dt, a = cfl_dt(q, dx, dy, gamma, CFL, R)
    it += 1

    if it % 50 == 0:
        print(f"{it:6d}  {t:10.4e}  {dt:10.4e}")

elapsed = time.perf_counter() - t_cpu
print(f"\nDone.  {it} iterations  |  CPU time: {elapsed:.2f} s")

# ============================================================
# Post-processing
# ============================================================
q_phys = q[R:ny-R, R:nx-R, :]
r  = q_phys[..., 0]
u  = q_phys[..., 1] / r
v  = q_phys[..., 2] / r
E  = q_phys[..., 3]
p  = (gamma - 1.0)*(E - 0.5*r*(u**2 + v**2))
c  = np.sqrt(np.maximum(gamma*p/r, 0.0))
U  = np.sqrt(u**2 + v**2)
M  = U / c
ss = np.log(p / r**gamma)          # dimensionless specific entropy
e  = p / ((gamma - 1.0)*r)         # specific internal energy

# ============================================================
# Final figure
# ============================================================
n_contours = 30
fig, axes = plt.subplots(2, 3, figsize=(14, 9))
titles  = [r'Density $\rho$', r'Speed $|\mathbf{u}|$', r'Pressure $p$',
           r'Entropy $\ln(p/\rho^\gamma)$', r'Mach number', r'Internal energy']
fields  = [r, U, p, ss, M, e]
for ax, fld, ttl in zip(axes.flat, fields, titles):
    cs = ax.contour(x, y, fld, n_contours)
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(ttl)

axes[0, 0].set_title(
    f'FV {RECON_MTH}-{FLUX_MTH} | Config {IC} | t={t:.3f}')
axes[0, 1].set_title(f'{NX}×{NY}, CFL={CFL}')
plt.tight_layout()

if SAVE_FIG:
    fig.savefig(FIG_PATH, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {FIG_PATH}")
plt.close(fig)

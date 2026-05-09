"""
FV-WENO solver for the Double Mach Reflection problem.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

Solves the classic Double Mach Reflection test problem of Woodward & Colella
(1984) using a Finite-Volume WENO scheme.  The domain is [0, 4] × [0, 1]
with a Mach-10 planar shock initially inclined at 30° at x = 1/6.

Boundary conditions
-------------------
- Left / right: zero-gradient (Neumann) outflow.
- Bottom:       slip wall (Neumann before the wedge, reflective after x = 1/6).
- Top:          time-dependent moving-shock Dirichlet / Neumann.

Usage
-----
    python run_fv_dmr.py

Edit the "Parameters" block to change resolution or scheme.

References
----------
[1] Woodward, P. and Colella, P., "The numerical simulation of two-dimensional
    fluid flow with strong shocks," J. Comput. Phys. 54, 115–173, 1984.
[2] Toro, E. F., "Riemann Solvers and Numerical Methods for Fluid Dynamics,"
    Springer, 2nd ed., 1999.
"""

import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from weno2d import fv_weno_rhs, ic_double_mach_reflection
from weno2d.reconstructions import stencil_radius
from weno2d.time_integration import ssprk3_step, cfl_dt

# ============================================================
# Parameters
# ============================================================
CFL       = 0.475        # CFL number (slightly conservative for DMR)
T_END     = 0.20         # final simulation time
NX        = 240          # physical cells in x
NY        = 60           # physical cells in y
N_DOF     = 5            # degrees of freedom (ideal diatomic gas)
FLUX_MTH  = 'HLLC'       # 'LF', 'LLF', 'ROE', 'HLLE', 'HLLC'
RECON_MTH = 'WENO5'      # 'WENO5', 'WENO7', 'Poly5', 'Poly7'
SAVE_FIG  = True
FIG_PATH  = 'dmr_result.png'

# ============================================================
# Derived parameters
# ============================================================
gamma = (N_DOF + 2.0) / N_DOF

# ============================================================
# Grid
# ============================================================
Lx = 4.0;  dx = Lx / NX
Ly = 1.0;  dy = Ly / NY
xc = np.linspace(dx/2.0, Lx - dx/2.0, NX)
yc = np.linspace(dy/2.0, Ly - dy/2.0, NY)
x, y = np.meshgrid(xc, yc)

# ============================================================
# Initial condition
# ============================================================
(r0, u0, v0, p0,
 preshock, postshock, shock_speed) = ic_double_mach_reflection(NX, NY, gamma)

E0 = p0/(gamma - 1.0) + 0.5*r0*(u0**2 + v0**2)
Q0 = np.stack([r0, r0*u0, r0*v0, E0], axis=-1)

# Wedge position in cell units (0-indexed into physical domain)
x_wedge        = 1.0/6.0
mesh_wedge_pos = x_wedge / dx      # used to locate the wall/wedge junction

# ============================================================
# Ghost-cell padding
# ============================================================
R  = stencil_radius(RECON_MTH)
nx = NX + 2*R;  ny = NY + 2*R
q  = np.zeros((ny, nx, 4))
q[R:ny-R, R:nx-R, :] = Q0

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
        test='DMR',
        gamma=gamma,
        mesh_wedge_position=mesh_wedge_pos,
        preshock=preshock,
        postshock=postshock,
        shock_speed=shock_speed,
    )

# ============================================================
# Time loop
# ============================================================
t     = 0.0
it    = 0
t_cpu = time.perf_counter()

print(f"\nFV-WENO DMR  |  {RECON_MTH}-{FLUX_MTH}  |  "
      f"{NX}×{NY} cells  |  t_end={T_END}")
print(f"{'iter':>6s}  {'t':>10s}  {'dt':>10s}")

while t < T_END:
    dt = min(dt, T_END - t)
    q  = ssprk3_step(q, dt, rhs)
    t += dt

    q_int = q[R:ny-R, R:nx-R, :]
    r_int = q_int[..., 0]
    p_int = (gamma - 1.0)*(q_int[..., 3]
             - 0.5*(q_int[..., 1]**2 + q_int[..., 2]**2) / r_int)
    if np.any(p_int < 0.0):
        raise RuntimeError(f"Negative pressure at t={t:.4e}.  "
                           "Reduce CFL or switch to a more dissipative flux.")

    dt, a = cfl_dt(q, dx, dy, gamma, CFL, R)
    it += 1
    if it % 20 == 0:
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
U  = np.sqrt(u**2 + v**2)
e  = p / ((gamma - 1.0)*r)

# ============================================================
# Final figures
# ============================================================
region = [0.0, 3.0, 0.0, 1.0]
n_contours = 30

fig, axes = plt.subplots(2, 2, figsize=(14, 7))
for ax, fld, ttl in zip(
        axes.flat,
        [r, U, p, e],
        [r'Density $\rho$', r'Speed $|\mathbf{u}|$',
         'Pressure $p$', 'Internal energy $e$']):
    cs = ax.contour(x, y, fld, n_contours)
    ax.axis(region)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(ttl)

axes[0, 0].set_title(f'FV {RECON_MTH}-{FLUX_MTH} | DMR | t={t:.3f}')
axes[0, 1].set_title(f'{NX}×{NY}, CFL={CFL}')
plt.tight_layout()

if SAVE_FIG:
    fig.savefig(FIG_PATH, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {FIG_PATH}")
plt.close(fig)

# Density contours from 1.731 to 20.92 (classical DMR benchmark range)
fig2, ax2 = plt.subplots(figsize=(14, 4))
v_levels = np.linspace(1.731, 20.92, n_contours)
ax2.contour(x, y, r, v_levels)
ax2.axis(region)
ax2.set_xlabel('x'); ax2.set_ylabel('y')
ax2.set_title(f'{n_contours} density contours 1.731 – 20.92 | '
              f'FV {RECON_MTH}-{FLUX_MTH} | t={t:.3f}')
plt.tight_layout()
if SAVE_FIG:
    fig2.savefig(FIG_PATH.replace('.png', '_density.png'), dpi=150, bbox_inches='tight')
plt.close(fig2)

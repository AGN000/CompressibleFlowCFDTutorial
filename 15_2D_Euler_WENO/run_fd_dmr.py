"""
FD-WENO solver for the Double Mach Reflection problem.

Modified from Manuel A. Diaz's WENO code (github.com/wme7/WENO).

Uses the Finite-Difference WENO method (Shu 2009) with Lax–Friedrichs or
local-Lax–Friedrichs flux splitting.

Usage
-----
    python run_fd_dmr.py

References
----------
[1] Woodward, P. and Colella, P., "The numerical simulation of two-dimensional
    fluid flow with strong shocks," J. Comput. Phys. 54, 115–173, 1984.
[2] Shu, C.-W., "High order weighted essentially non-oscillatory schemes for
    convection dominated problems," SIAM Review 51(1), 82–126, 2009.
"""

import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from weno2d import fd_weno_rhs, ic_double_mach_reflection
from weno2d.reconstructions import stencil_radius
from weno2d.time_integration import ssprk3_step, cfl_dt

# ============================================================
# Parameters
# ============================================================
CFL        = 0.475
T_END      = 0.20
NX         = 241
NY         = 61
N_DOF      = 5
SPLIT_MTH  = 'LF'         # 'LF' or 'LLF'
RECON_MTH  = 'WENO5'      # 'WENO5', 'WENO7', 'Poly5', 'Poly7'
SAVE_FIG   = True
FIG_PATH   = 'fd_dmr_result.png'

# ============================================================
# Derived parameters
# ============================================================
gamma = (N_DOF + 2.0) / N_DOF

# ============================================================
# Grid
# ============================================================
Lx = 4.0;  dx = Lx / (NX - 1)
Ly = 1.0;  dy = Ly / (NY - 1)
xc = np.linspace(0.0, Lx, NX)
yc = np.linspace(0.0, Ly, NY)
x, y = np.meshgrid(xc, yc)

# ============================================================
# Initial condition
# ============================================================
(r0, u0, v0, p0,
 preshock, postshock, shock_speed) = ic_double_mach_reflection(NX, NY, gamma)

E0 = p0/(gamma - 1.0) + 0.5*r0*(u0**2 + v0**2)
Q0 = np.stack([r0, r0*u0, r0*v0, E0], axis=-1)

x_wedge        = 1.0/6.0
mesh_wedge_pos = x_wedge / dx

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
# Bind RHS
# ============================================================
def rhs(q_in):
    return fd_weno_rhs(
        q_in, a, nx, ny, dx, dy, t,
        split_method=SPLIT_MTH,
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

print(f"\nFD-WENO DMR  |  {RECON_MTH}-{SPLIT_MTH}  |  "
      f"{NX}×{NY} cells  |  t_end={T_END}")
print(f"{'iter':>6s}  {'t':>10s}  {'dt':>10s}")

while t < T_END:
    dt = min(dt, T_END - t)
    q  = ssprk3_step(q, dt, rhs)
    t += dt
    dt, a = cfl_dt(q, dx, dy, gamma, CFL, R)
    it += 1
    if it % 20 == 0:
        print(f"{it:6d}  {t:10.4e}  {dt:10.4e}")

elapsed = time.perf_counter() - t_cpu
print(f"\nDone.  {it} iterations  |  CPU time: {elapsed:.2f} s")

# ============================================================
# Post-processing and figure
# ============================================================
q_phys = q[R:ny-R, R:nx-R, :]
r  = q_phys[..., 0]
u  = q_phys[..., 1] / r
v  = q_phys[..., 2] / r
E  = q_phys[..., 3]
p  = (gamma - 1.0)*(E - 0.5*r*(u**2 + v**2))
U  = np.sqrt(u**2 + v**2)

region = [0.0, 3.0, 0.0, 1.0]
n_contours = 30

fig, axes = plt.subplots(2, 2, figsize=(14, 7))
for ax, fld, ttl in zip(axes.flat, [r, U, p, p/((gamma-1)*r)],
                         [r'$\rho$', r'$|\mathbf{u}|$', '$p$', '$e$']):
    ax.contour(x, y, fld, n_contours)
    ax.axis(region); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_title(ttl)

axes[0, 0].set_title(f'FD {RECON_MTH}-{SPLIT_MTH} | DMR | t={t:.3f}')
axes[0, 1].set_title(f'{NX}×{NY}, CFL={CFL}')
plt.tight_layout()
if SAVE_FIG:
    fig.savefig(FIG_PATH, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {FIG_PATH}")
plt.close(fig)

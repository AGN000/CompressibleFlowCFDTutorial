# 15 — 2D Compressible Euler Equations: WENO Solver

A fully vectorised Python solver for the **2D compressible Euler equations** using
Weighted Essentially Non-Oscillatory (WENO) reconstruction with both Finite-Volume
(FV) and Finite-Difference (FD) formulations.

Modified from [Manuel A. Diaz's WENO code](https://github.com/wme7/WENO).

---

## Governing Equations

$$\frac{\partial \mathbf{q}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} = 0$$

with conserved state $\mathbf{q} = [\rho,\, \rho u,\, \rho v,\, E]^\top$ and ideal-gas
closure $p = (\gamma - 1)\!\left(E - \tfrac{1}{2}\rho(u^2 + v^2)\right)$.

---

## Package structure

```
weno2d/
├── reconstructions.py      # WENO5 / WENO7 / Poly5 / Poly7 (x and y directions)
├── flux_solvers.py         # LF, LLF/Rusanov, Roe, HLLE, HLLC  (FV)
│                           # LF / LLF flux splitting             (FD)
├── boundary_conditions.py  # Periodic, Neumann outflow, DMR moving-shock BC
├── rhs_fv.py               # Finite-volume WENO residual
├── rhs_fd.py               # Finite-difference WENO residual
├── initial_conditions.py   # DMR IC (Rankine–Hugoniot) + 19 Riemann configs
└── time_integration.py     # SSP-RK3 (Shu–Osher) + CFL time-step control
```

---

## Driver scripts

| Script | Problem | Method |
|---|---|---|
| `run_fv_riemann2d.py` | 2D Riemann problem (any of 19 configs) | FV-WENO |
| `run_fv_dmr.py` | Double Mach Reflection | FV-WENO |
| `run_fd_dmr.py` | Double Mach Reflection | FD-WENO |

Edit the **Parameters** block at the top of each script to change resolution,
scheme, flux method, or configuration.

### Quick start

```bash
pip install numpy matplotlib
python run_fv_riemann2d.py   # Config 6, 200×200, HLLC, t=0.3
python run_fv_dmr.py         # DMR 240×60, HLLC, t=0.2
python run_fd_dmr.py         # DMR 241×61, LF,   t=0.2
python test_weno2d.py        # 6-test validation suite
```

---

## Numerical methods

### Reconstruction
- **WENO5** — 5th-order, stencil radius R=3 (default)
- **WENO7** — 7th-order, stencil radius R=4
- **Poly5 / Poly7** — linear (no nonlinear weights), same stencil

### Riemann solvers (FV only)
| Name | Description |
|---|---|
| `LF` | Global Lax–Friedrichs |
| `LLF` | Local Lax–Friedrichs / Rusanov |
| `ROE` | Roe with Harten entropy fix |
| `HLLE` | Einfeldt wave-speed estimate |
| `HLLC` | Two-contact HLLC (Toro) |

### Flux splitting (FD only)
- `LF` — global Lax–Friedrichs $F^\pm = (F \pm a_{\max} q)/2$
- `LLF` — local Lax–Friedrichs (cell-by-cell spectral radius)

### Time integration
SSP-RK3 (Shu–Osher convex-combination form):

$$q^{(1)} = q^n + \Delta t\, L(q^n)$$
$$q^{(2)} = \tfrac{3}{4}q^n + \tfrac{1}{4}(q^{(1)} + \Delta t\, L(q^{(1)}))$$
$$q^{n+1} = \tfrac{1}{3}q^n + \tfrac{2}{3}(q^{(2)} + \Delta t\, L(q^{(2)}))$$

---

## Test cases

### Double Mach Reflection (Woodward & Colella 1984)
- Mach-10 planar shock inclined at 30° on domain $[0,4]\times[0,1]$
- Moving-shock Dirichlet BC at the top, reflective wall at the bottom
- Benchmark: 30 density contours between 1.731 and 20.92

### 2D Riemann problems (Kurganov & Tadmor 2002)
- 19 standard configurations on the unit square $[0,1]^2$
- Initial data: four constant states in the four quadrants

---

## Validation

```
python test_weno2d.py
```

Six automated checks:

1. WENO5 reconstruction symmetry
2. 5th-order convergence on smooth data (orders ≥ 4.5)
3. Rusanov flux consistency $F(q, q) = F_{\text{phys}}(q)$
4. DMR IC satisfies Rankine–Hugoniot conditions
5. Constant flow produces zero FV residual
6. One SSP-RK3 step stays finite with positive density and pressure

---

## Dependencies

- Python ≥ 3.9
- NumPy
- Matplotlib (figures only)

---

## References

1. Woodward, P. and Colella, P., *The numerical simulation of two-dimensional fluid
   flow with strong shocks*, J. Comput. Phys. **54**, 115–173, 1984.
2. Shu, C.-W., *High order weighted essentially non-oscillatory schemes for
   convection dominated problems*, SIAM Review **51**(1), 82–126, 2009.
3. Kurganov, A. and Tadmor, E., *Solution of two-dimensional Riemann problems for
   gas dynamics without Riemann problem solvers*, Numer. Methods PDE **18**(5),
   584–608, 2002.
4. Toro, E. F., *Riemann Solvers and Numerical Methods for Fluid Dynamics*,
   Springer, 2nd ed., 1999.

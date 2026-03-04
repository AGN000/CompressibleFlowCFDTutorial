# 🚀 Compressible Flow CFD Tutorial (Python)

A structured, research-oriented tutorial series on **compressible flow CFD**, built from first principles using **pure Python**.

This repository builds intuition first — starting from linear convection and progressing toward high-speed compressible Euler equations, shock capturing, Riemann solvers, and high-order schemes.

---

## 🎯 Objectives

This tutorial series aims to:

- Build strong intuition for conservation laws
- Understand instability of central schemes
- Introduce the Finite Volume Method (FVM)
- Implement shock-capturing methods
- Explore higher-order schemes (TVD, WENO)
- Extend from 1D to 2D compressible Euler equations
- Provide a foundation for hybrid CFD–ML research

---

# 📘 Tutorial Roadmap

---

## 1️⃣ Linear Convection Equation (Foundation)

We begin with the 1D linear convection equation:

$$
\frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} = 0
$$

with periodic boundary conditions.

### Topics Covered

- FTCS, FTBS, FTFS schemes
- CFL condition and stability
- Numerical diffusion vs dispersion
- Why central schemes fail
- Need for upwinding
- Error comparison with exact solution
- Python animations

---

## 2️⃣ Finite Volume Method (FVM)

We transition to conservation form:

$$
\frac{\partial}{\partial t} \int_{V_i} U \, dV 
+ \int_{\partial V_i} F(U)\cdot n \, dS = 0
$$

### Topics

- Integral form of conservation laws
- Cell-averaged variables
- Numerical flux functions
- First-Order Upwind (FOU)
- Conservative discretization
- Clean Python implementation

---

## 3️⃣ 1D Euler Equations

$$
\frac{\partial}{\partial t}
\begin{bmatrix}
\rho \\
\rho u \\
E
\end{bmatrix}
+
\frac{\partial}{\partial x}
\begin{bmatrix}
\rho u \\
\rho u^2 + p \\
u(E+p)
\end{bmatrix}
= 0
$$

### Topics

- Conservative variables
- Primitive ↔ Conservative transformation
- Sod shock tube problem
- Numerical flux construction

---

## 4️⃣ Riemann Solvers

- Exact Riemann Solver
- Roe Solver
- HLL
- HLLC

Wave structures:
- Shock
- Rarefaction
- Contact discontinuity

---

## 5️⃣ Higher-Order Shock Capturing

- MUSCL reconstruction
- TVD limiters (Minmod, Van Leer, Superbee)
- Total Variation concept
- WENO schemes
- Non-oscillatory high-order reconstruction

---

## 6️⃣ Extension to 2D Compressible Flows

- 2D Euler equations
- Dimensional splitting
- 2D Riemann solvers
- Structured grid implementation

---

# 🧠 Philosophy

Many CFD resources:

- Jump directly to complex C++ solvers
- Or skip physical intuition

This repository emphasizes:

✔ Flux interpretation  
✔ Stability understanding  
✔ Clean Python implementations  
✔ Visualization and diagnostics  

Intuition → Stability → Conservation → Shock Capturing → High Order

---

# 🔬 Research Motivation

This tutorial series also builds a foundation for:

- Adaptive viscosity modeling
- Physics-Informed Neural Networks (PINNs)
- Differentiable CFD solvers
- Hybrid ML + FVM frameworks

The long-term goal:

**Bridge classical numerical methods with physics-informed learning.**

---


# ============================================================
# Neural Artificial Viscosity Learning
# ------------------------------------------------------------
# Learn ν(x,t) directly from local solution stencil
#
# PDE:
#
#     u_t + c u_x = 0
#
# Numerical Scheme:
#
#     u^{n+1}
#     =
#     u^n
#     -
#     c dt/(2dx) (u_{i+1}-u_{i-1})
#     +
#     ν_i dt/dx² (u_{i+1}-2u_i+u_{i-1})
#
# Neural Network learns ν_i.
#
# Loss:
#
#     MSE(predicted trajectory, exact trajectory)
#
# Optimization:
#
#     Adam + L-BFGS
#
# ============================================================

import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt

# ============================================================
# Device
# ============================================================

device = torch.device(
    "cuda" if torch.cuda.is_available() else "cpu"
)

print("Device:", device)

# ============================================================
# Problem parameters
# ============================================================

L = 1.0

N = 100

dx = L / N

c = 1.0

dt = 0.001

total_steps = 150

adam_epochs = 1000

# Allow small negative viscosity

nu_min = -0.005
nu_max = 0.02

# ============================================================
# Grid
# ============================================================

x = np.linspace(
    0.0,
    L,
    N,
    endpoint=False
)

# ============================================================
# Initial condition
# Hat function
# ============================================================

u0_np = np.zeros(
    N,
    dtype=np.float32
)

u0_np[
    (x > 0.4) &
    (x < 0.6)
] = 1.0

u0 = torch.tensor(
    u0_np,
    dtype=torch.float32,
    device=device
)

# ============================================================
# Exact solution trajectory
# ============================================================

exact_traj = []

for n in range(1, total_steps + 1):

    shift = int(
        c * n * dt / dx
    )

    u_exact = np.roll(
        u0_np,
        shift
    )

    exact_traj.append(
        u_exact
    )

exact_traj = np.array(
    exact_traj,
    dtype=np.float32
)

exact_traj_torch = torch.tensor(
    exact_traj,
    dtype=torch.float32,
    device=device
)

# ============================================================
# Neural viscosity model
# ============================================================

class NuNet(nn.Module):

    """
    Input:
        [u_{i-1}, u_i, u_{i+1}]

    Output:
        ν_i
    """

    def __init__(self):

        super().__init__()

        self.net = nn.Sequential(

            nn.Linear(3, 64),
            nn.Tanh(),

            nn.Linear(64, 64),
            nn.Tanh(),

            nn.Linear(64, 64),
            nn.Tanh(),

            nn.Linear(64, 1)

        )

    def forward(self, um, u, up):

        inp = torch.stack(
            [um, u, up],
            dim=1
        )

        nu = self.net(inp)

        # bounded output

        nu = 0.02 * torch.tanh(nu)

        return nu


net = NuNet().to(device)

# ============================================================
# Differentiable rollout
# ============================================================

def rollout(model):

    u = u0.clone()

    trajectory = []

    nu_history = []

    for _ in range(total_steps):

        up = torch.roll(
            u,
            shifts=-1,
            dims=0
        )

        um = torch.roll(
            u,
            shifts=1,
            dims=0
        )

        # ---------------------------------
        # Neural viscosity
        # ---------------------------------

        nu = model(
            um,
            u,
            up
        ).squeeze()

        nu = torch.clamp(
            nu,
            nu_min,
            nu_max
        )

        # ---------------------------------
        # Central advection
        # ---------------------------------

        adv = -(
            c * dt /
            (2.0 * dx)
        ) * (up - um)

        # ---------------------------------
        # Diffusion term
        # ---------------------------------

        lap = (
            up
            - 2.0 * u
            + um
        )

        diff = (
            nu
            * dt
            / dx**2
        ) * lap

        # ---------------------------------
        # Update
        # ---------------------------------

        u = u + adv + diff

        trajectory.append(u)

        nu_history.append(nu)

    trajectory = torch.stack(
        trajectory,
        dim=0
    )

    nu_history = torch.stack(
        nu_history,
        dim=0
    )

    return trajectory, nu_history

# ============================================================
# Loss function
# ============================================================

def compute_loss(model):

    pred_traj, nu_traj = rollout(model)

    solution_loss = torch.mean(
        (
            pred_traj
            -
            exact_traj_torch
        )**2
    )

    # Regularization

    nu_reg = 1e-4 * torch.mean(
        nu_traj**2
    )

    # Penalize large negative ν

    neg_reg = 1e-3 * torch.mean(
        torch.relu(
            -nu_traj
        )**2
    )

    loss = (
        solution_loss
        + nu_reg
        + neg_reg
    )

    return (
        loss,
        solution_loss,
        pred_traj,
        nu_traj
    )

# ============================================================
# Adam optimization
# ============================================================

optimizer = torch.optim.Adam(
    net.parameters(),
    lr=1e-3
)

print("\nStarting Adam...\n")

for epoch in range(adam_epochs):

    optimizer.zero_grad()

    loss, sol_loss, _, _ = compute_loss(
        net
    )

    loss.backward()

    torch.nn.utils.clip_grad_norm_(
        net.parameters(),
        1.0
    )

    optimizer.step()

    if epoch % 100 == 0:

        print(
            f"Epoch {epoch:5d} | "
            f"Total={loss.item():.6e} | "
            f"Solution={sol_loss.item():.6e}"
        )

# ============================================================
# L-BFGS refinement
# ============================================================

print("\nStarting L-BFGS...\n")

lbfgs = torch.optim.LBFGS(
    net.parameters(),
    lr=1.0,
    max_iter=50,
    max_eval=500,
    history_size=100,
    line_search_fn="strong_wolfe"
)

def closure():

    lbfgs.zero_grad()

    loss, _, _, _ = compute_loss(
        net
    )

    loss.backward()

    return loss

lbfgs.step(
    closure
)

# ============================================================
# Final evaluation
# ============================================================

final_loss, final_sol_loss, pred_traj, nu_traj = compute_loss(
    net
)

print("\nFinal loss =", final_loss.item())
print(
    "Final solution loss =",
    final_sol_loss.item()
)

# ============================================================
# Convert to numpy
# ============================================================

pred_traj_np = (
    pred_traj
    .detach()
    .cpu()
    .numpy()
)

nu_traj_np = (
    nu_traj
    .detach()
    .cpu()
    .numpy()
)

u_final_pred = pred_traj_np[-1]

u_final_exact = exact_traj[-1]

# ============================================================
# L2 error
# ============================================================

l2_error = np.sqrt(
    np.mean(
        (
            u_final_pred
            -
            u_final_exact
        )**2
    )
)

print(
    "\nFinal L2 Error =",
    l2_error
)

# ============================================================
# Plot solution
# ============================================================

plt.figure(
    figsize=(10,5)
)

plt.plot(
    x,
    u_final_exact,
    'k--',
    lw=3,
    label='Exact'
)

plt.plot(
    x,
    u_final_pred,
    'r',
    lw=2,
    label='NN Learned ν'
)

plt.xlabel("x")
plt.ylabel("u")

plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# ============================================================
# Plot ν(x,t)
# ============================================================

MU = nu_traj_np.T

times = (
    np.arange(
        1,
        total_steps + 1
    ) * dt
)

T, X = np.meshgrid(
    times,
    x
)

plt.figure(
    figsize=(12,5)
)

plt.contourf(
    X,
    T,
    MU,
    levels=100,
    cmap="magma"
)

plt.colorbar(
    label="ν(x,t)"
)

plt.xlabel("x")
plt.ylabel("t")

plt.tight_layout()
plt.show()

# ============================================================
# Error contour
# ============================================================

ERR = (
    pred_traj_np
    -
    exact_traj
)

plt.figure(
    figsize=(12,5)
)

plt.contourf(
    X,
    T,
    ERR.T,
    levels=100
)

plt.colorbar(
    label="u - u_exact"
)

plt.xlabel("x")
plt.ylabel("t")

plt.tight_layout()
plt.show()

print("\nDone.")

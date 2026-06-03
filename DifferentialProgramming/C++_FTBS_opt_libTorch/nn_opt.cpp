// ============================================================
// Neural Artificial Viscosity Learning - C++ (LibTorch)
// ------------------------------------------------------------
// Learn nu(x,t) directly from local solution stencil
//
// PDE:
//     u_t + c u_x = 0
//
// Numerical Scheme:
//     u^{n+1} = u^n - c*dt/(2*dx)*(u_{i+1}-u_{i-1})
//               + nu_i*dt/dx^2*(u_{i+1}-2*u_i+u_{i-1})
//
// Neural Network learns nu_i.
//
// Loss: MSE(predicted trajectory, exact trajectory)
// Optimization: Adam + L-BFGS
// ============================================================

#include <torch/torch.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

// ============================================================
// Device
// ============================================================

auto device = torch::cuda::is_available()
    ? torch::kCUDA
    : torch::kCPU;

// ============================================================
// Problem parameters
// ============================================================

const double L = 1.0;
const int N = 100;
const double dx = L / N;
const double c = 1.0;
const double dt = 0.001;
const int total_steps = 150;
const int adam_epochs = 1000;

const double nu_min = -0.005;
const double nu_max = 0.02;

// ============================================================
// Neural viscosity model
// ============================================================

struct NuNet : torch::nn::Module {
    torch::nn::Sequential net{nullptr};

    NuNet() {
        net = register_module("net", torch::nn::Sequential(
            torch::nn::Linear(3, 64),
            torch::nn::Tanh(),
            torch::nn::Linear(64, 64),
            torch::nn::Tanh(),
            torch::nn::Linear(64, 64),
            torch::nn::Tanh(),
            torch::nn::Linear(64, 1)
        ));
    }

    torch::Tensor forward(torch::Tensor um, torch::Tensor u, torch::Tensor up) {
        auto inp = torch::stack({um, u, up}, 1);
        auto nu = net->forward(inp);
        nu = 0.02 * torch::tanh(nu);
        return nu;
    }
};

// ============================================================
// Differentiable rollout
// ============================================================

struct RolloutResult {
    torch::Tensor trajectory;
    torch::Tensor nu_history;
};

RolloutResult rollout(
    NuNet& model,
    const torch::Tensor& u0,
    const torch::Tensor& exact_traj_torch
) {
    auto u = u0.clone();

    std::vector<torch::Tensor> trajectory_vec;
    std::vector<torch::Tensor> nu_history_vec;
    trajectory_vec.reserve(total_steps);
    nu_history_vec.reserve(total_steps);

    for (int i = 0; i < total_steps; ++i) {
        auto up = torch::roll(u, -1, 0);
        auto um = torch::roll(u, 1, 0);

        auto nu = model.forward(um, u, up).squeeze();
        nu = torch::clamp(nu, nu_min, nu_max);

        auto adv = -(c * dt / (2.0 * dx)) * (up - um);

        auto lap = up - 2.0 * u + um;
        auto diff = (nu * dt / (dx * dx)) * lap;

        u = u + adv + diff;

        trajectory_vec.push_back(u);
        nu_history_vec.push_back(nu);
    }

    auto trajectory = torch::stack(trajectory_vec, 0);
    auto nu_history = torch::stack(nu_history_vec, 0);

    return {trajectory, nu_history};
}

// ============================================================
// Loss function
// ============================================================

struct LossResult {
    torch::Tensor loss;
    torch::Tensor solution_loss;
    torch::Tensor trajectory;
    torch::Tensor nu_history;
};

LossResult compute_loss(
    NuNet& model,
    const torch::Tensor& u0,
    const torch::Tensor& exact_traj_torch
) {
    auto [pred_traj, nu_traj] = rollout(model, u0, exact_traj_torch);

    auto solution_loss = torch::mean(torch::pow(pred_traj - exact_traj_torch, 2));

    auto nu_reg = 1e-4 * torch::mean(torch::pow(nu_traj, 2));

    auto neg_reg = 1e-3 * torch::mean(torch::pow(torch::relu(-nu_traj), 2));

    auto loss = solution_loss + nu_reg + neg_reg;

    return {loss, solution_loss, pred_traj, nu_traj};
}

// ============================================================
// Write 2D field to CSV
// ============================================================

void write_csv(
    const std::string& filename,
    const torch::Tensor& data,
    const torch::Tensor& x_grid,
    const torch::Tensor& t_grid
) {
    std::ofstream f(filename);
    f << "x,t,value\n";
    auto d = data.contiguous();
    auto xg = x_grid.contiguous();
    auto tg = t_grid.contiguous();
    auto data_a = d.accessor<float, 2>();
    auto x_a = xg.accessor<float, 2>();
    auto t_a = tg.accessor<float, 2>();
    for (int i = 0; i < data.size(0); ++i) {
        for (int j = 0; j < data.size(1); ++j) {
            f << x_a[i][j] << "," << t_a[i][j] << "," << data_a[i][j] << "\n";
        }
    }
}

int main() {
    std::cout << "Device: " << (device == torch::kCUDA ? "cuda" : "cpu") << std::endl;

    // ============================================================
    // Initial condition - hat function
    // ============================================================

    auto x = torch::linspace(0.0, L, N, torch::kFloat32);

    auto u0_np = torch::zeros({N}, torch::kFloat32);
    auto mask = (x > 0.4) & (x < 0.6);
    u0_np.masked_fill_(mask, 1.0f);

    auto u0 = u0_np.to(device);

    // ============================================================
    // Exact solution trajectory
    // ============================================================

    std::vector<torch::Tensor> exact_vec;
    for (int n = 1; n <= total_steps; ++n) {
        int shift = static_cast<int>(c * n * dt / dx);
        exact_vec.push_back(torch::roll(u0_np, shift));
    }
    auto exact_traj = torch::stack(exact_vec, 0);
    auto exact_traj_torch = exact_traj.to(device);

    // ============================================================
    // Model
    // ============================================================

    NuNet net;
    net.to(device);

    // ============================================================
    // Adam optimization
    // ============================================================

    torch::optim::Adam adam_opt(net.parameters(), 1e-3);

    std::cout << "\nStarting Adam...\n" << std::endl;

    for (int epoch = 0; epoch < adam_epochs; ++epoch) {
        adam_opt.zero_grad();

        auto [loss, sol_loss, _, __] = compute_loss(net, u0, exact_traj_torch);

        loss.backward();

        torch::nn::utils::clip_grad_norm_(net.parameters(), 1.0);

        adam_opt.step();

        if (epoch % 100 == 0) {
            std::cout << "Epoch " << std::setw(5) << epoch << " | "
                      << "Total=" << std::scientific << loss.item<double>() << " | "
                      << "Solution=" << std::scientific << sol_loss.item<double>() << std::endl;
        }
    }

    // ============================================================
    // L-BFGS refinement
    // ============================================================

    std::cout << "\nStarting L-BFGS...\n" << std::endl;

    torch::optim::LBFGS lbfgs_opt(net.parameters(),
        torch::optim::LBFGSOptions()
            .lr(1.0)
            .max_iter(5)
            .max_eval(50)
            .history_size(10)
            .line_search_fn("strong_wolfe"));

    auto closure = [&]() -> torch::Tensor {
        lbfgs_opt.zero_grad();
        auto [loss, _, __, ___] = compute_loss(net, u0, exact_traj_torch);
        loss.backward();
        return loss;
    };

    lbfgs_opt.step(closure);

    // ============================================================
    // Final evaluation
    // ============================================================

    auto [final_loss, final_sol_loss, pred_traj, nu_traj] =
        compute_loss(net, u0, exact_traj_torch);

    std::cout << "\nFinal loss = " << final_loss.item<double>() << std::endl;
    std::cout << "Final solution loss = " << final_sol_loss.item<double>() << std::endl;

    // ============================================================
    // Convert to CPU / numpy equivalents
    // ============================================================

    auto pred_traj_cpu = pred_traj.cpu();
    auto nu_traj_cpu = nu_traj.cpu();
    auto exact_traj_cpu = exact_traj.cpu();

    auto u_final_pred = pred_traj_cpu[total_steps - 1];
    auto u_final_exact = exact_traj_cpu[total_steps - 1];

    // ============================================================
    // L2 error
    // ============================================================

    auto l2_error = torch::sqrt(torch::mean(torch::pow(u_final_pred - u_final_exact, 2)));
    std::cout << "\nFinal L2 Error = " << l2_error.item<double>() << std::endl;

    // ============================================================
    // Write results to CSV files for external plotting
    // ============================================================

    // Final time slice comparison
    {
        std::ofstream f("final_comparison.csv");
        f << "x,exact,predicted\n";
        auto x_a = x.accessor<float, 1>();
        auto exact_a = u_final_exact.accessor<float, 1>();
        auto pred_a = u_final_pred.accessor<float, 1>();
        for (int i = 0; i < N; ++i) {
            f << x_a[i] << "," << exact_a[i] << "," << pred_a[i] << "\n";
        }
    }

    // nu(x,t) field
    {
        auto times = torch::arange(1, total_steps + 1, torch::kFloat32) * dt;
        auto T = times.unsqueeze(1).expand({total_steps, N});
        auto X = x.unsqueeze(0).expand({total_steps, N});
        auto MU = nu_traj_cpu.transpose(0, 1);

        write_csv("nu_field.csv", MU, X.transpose(0, 1), T.transpose(0, 1));
    }

    // Error field
    {
        auto ERR = pred_traj_cpu - exact_traj_cpu;
        auto times = torch::arange(1, total_steps + 1, torch::kFloat32) * dt;
        auto T = times.unsqueeze(1).expand({total_steps, N});
        auto X = x.unsqueeze(0).expand({total_steps, N});

        write_csv("error_field.csv", ERR, X.transpose(0, 1), T.transpose(0, 1));
    }

    std::cout << "\nOutput written to: final_comparison.csv, nu_field.csv, error_field.csv" << std::endl;
    std::cout << "\nDone." << std::endl;

    return 0;
}

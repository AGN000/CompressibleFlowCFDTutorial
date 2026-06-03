# Neural Artificial Viscosity Learning — C++

Learn `ν(x,t)` directly from a local solution stencil for the 1D advection PDE `u_t + c u_x = 0`.

## Requirements

- **C++ compiler** with C++17 support (GCC ≥ 9, Clang ≥ 10)
- **CMake** ≥ 3.18 (optional — a direct `g++` build is also supported)
- **LibTorch / PyTorch** ≥ 2.0 (C++ API, GPU or CPU)

## Installing LibTorch

### Option A — Use an existing PyTorch installation (recommended)

If you already have PyTorch installed via `pip` or `conda`, it ships with the C++ headers and libraries.

Find the cmake prefix path:

```bash
python3 -c "import torch; print(torch.utils.cmake_prefix_path)"
```

Then configure CMake with that path:

```bash
cmake -B build \
  -DCMAKE_PREFIX_PATH="$(python3 -c 'import torch; print(torch.utils.cmake_prefix_path)')" \
  -DCMAKE_DISABLE_FIND_PACKAGE_CUDA=ON   # only if nvcc causes issues with your compiler
```

### Option B — Download standalone LibTorch

1. Go to https://pytorch.org/get-started/locally/
2. Select your platform and choose **C++ / Java** (not pip/conda)
3. Copy the download link for the **cxx11 ABI** version (compatible with most GCC installations)
4. Download and extract:

```bash
wget <download-url> -O libtorch.zip
unzip libtorch.zip
```

Then configure with:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/libtorch
```

## Building

### With CMake

```bash
cmake -B build -DCMAKE_PREFIX_PATH=<path-to-libtorch>
cmake --build build
```

### Without CMake (direct g++)

For a conda/miniconda environment with PyTorch installed:

```bash
TORCH_DIR=$(python3 -c "import torch; import os; print(os.path.dirname(torch.__file__))")
CUDA_DIR=$CONDA_PREFIX/targets/x86_64-linux
CONDA_LIB=$CONDA_PREFIX/lib

g++ -std=c++17 \
  -I$TORCH_DIR/include \
  -I$TORCH_DIR/include/torch/csrc/api/include \
  -I$CUDA_DIR/include \
  -L$TORCH_DIR/lib \
  -L$CUDA_DIR/lib \
  -L$CONDA_LIB \
  -Wl,-rpath,$TORCH_DIR/lib \
  -Wl,-rpath,$CUDA_DIR/lib \
  -Wl,-rpath,$CONDA_LIB \
  nn_opt.cpp \
  -ltorch -ltorch_cpu -ltorch_cuda -lc10 -lc10_cuda -lcudart \
  -lsleef -lprotobuf \
  -o nn_opt
```

## Running

```bash
./nn_opt
```

The program:
1. Runs **Adam** for 1000 epochs
2. Runs **L-BFGS** for 5 refinement iterations
3. Writes three CSV files: `final_comparison.csv`, `nu_field.csv`, `error_field.csv`
4. Prints the final loss and L2 error

## Visualising results

Install Python with `pandas` and `matplotlib`, then:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Final time slice
fc = pd.read_csv('final_comparison.csv')
plt.plot(fc['x'], fc['exact'], 'k--', label='Exact')
plt.plot(fc['x'], fc['predicted'], 'r', label='NN Learned nu')
plt.legend()
plt.savefig('final_comparison.png')

# nu field contour
nu = pd.read_csv('nu_field.csv')
X = nu['x'].values.reshape(100, 150)
T = nu['t'].values.reshape(100, 150)
MU = nu['value'].values.reshape(100, 150)
plt.contourf(X, T, MU, levels=100, cmap='magma')
plt.colorbar(label='nu(x,t)')
plt.savefig('nu_field.png')

# Error contour
err = pd.read_csv('error_field.csv')
ERR = err['value'].values.reshape(100, 150)
plt.contourf(X, T, ERR, levels=100)
plt.colorbar(label='u - u_exact')
plt.savefig('error_field.png')
```

## Parameters

Edit `nn_opt.cpp` to change:

| Parameter | Default | Description |
|---|---|---|
| `N` | 100 | Grid points |
| `c` | 1.0 | Advection speed |
| `dt` | 0.001 | Time step |
| `total_steps` | 150 | Rollout length |
| `adam_epochs` | 1000 | Adam iterations |
| `nu_min` / `nu_max` | -0.005 / 0.02 | Viscosity bounds |

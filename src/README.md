# Source Code - MATLAB Implementation

This folder contains all MATLAB source files implementing ADMM algorithms for inverse problems.

## File Organization

### Main Entry Point

- **`A_MAIN_2D.m`** (112 KB) — Master script coordinating all experiments
  - Loads test images
  - Generates degraded versions
  - Runs multiple ADMM variants
  - Computes and visualizes results
  - **Start here to run experiments**

### Data Generation

- **`B1_REST_GENERATE_DATA_2D.m`** — Generate restoration (deblurring) problem data
  - Creates blurred images with noise
  - Implements blur kernels

- **`B2_INPT_GENERATE_DATA_2D.m`** — Generate inpainting problem data
  - Creates masked images
  - Implements super-resolution masking

### Restoration Algorithms (Deblurring)

These solve the **restoration** problem: recover a blurred image with unknown noise.

| File | Regularization | Penalty | Method |
|------|---|---|---|
| `REST_TV_L1_U_ADMM.m` | Total Variation | L1 | Standard ADMM |
| `REST_TV_L2_DC_ADMM.m` | Total Variation | L2 | Dual Compositional |
| `REST_TIK_L1_U_ADMM.m` | Tikhonov | L1 | Standard ADMM |
| `REST_TIK_L2_U_DIR.m` | Tikhonov | L2 | Direct minimization |

### Inpainting Algorithms

These solve the **inpainting** problem: recover missing pixels (equivalent to super-resolution with sr-factor=2).

#### Standard Implementations

| File | Regularization | Penalty |
|------|---|---|
| `INPT_TV_L1_U_ADMM.m` | Total Variation | L1 |
| `INPT_TV_L2_U_ADMM.m` | Total Variation | L2 |
| `INPT_TIK_L1_U_ADMM.m` | Tikhonov | L1 |
| `INPT_TIK_L2_U_DIR.m` | Tikhonov | L2 |

#### L1.5 Regularization (Enhanced Formulations)

These implement sparsity penalties with exponent 1.5 — between L1 (sparse) and L2 (smooth), often providing better edge preservation.

| File | Regularization | Penalty | Method |
|------|---|---|---|
| `INPT_TV_L15_U_ADMM.m` | Total Variation | L1.5 | Standard ADMM |
| `INPT_TV_L15_U_ADMM_FAST.m` | Total Variation | L1.5 | **Fast ADMM** |
| `INPT_TIK_L15_U_ADMM.m` | Tikhonov | L1.5 | Standard ADMM |
| `INPT_TIK_L15_U_ADMM_FAST.m` | Tikhonov | L1.5 | **Fast ADMM** |

### Utility Functions

#### Signal Processing

- **`MASK_IMAGE_COL.m`** — Apply masking to images (for inpainting setup)
- **`compute_snr.m`** — Compute Signal-to-Noise Ratio
- **`compute_normalized_histogram.m`** — Compute normalized image histogram

#### Noise Models (PDFs)

These compute probability density functions for various noise distributions:

- **`gauss_pdf_1D.m`** — Gaussian distribution
- **`laplace_pdf_1D.m`** — Laplace distribution
- **`generalized_gauss_pdf_1D.m`** — Generalized Gaussian: p(x|α,β,μ) = (β/2αΓ(1/β)) exp(-(|x-μ|/α)^β)
- **`uniform_pdf_1D.m`** — Uniform distribution
- **`generate_WGG_realization.m`** — Generate Multiplicative White Gaussian Gaussian (MWGG) realizations

### Analysis Scripts

- **`LAB1_POINT1.m`** — Noise inference and characterization
  - Analyzes image sequence to estimate noise statistics
  - Fits distribution models (Gaussian, Laplace, Generalized Gaussian, MWGG)
  - Generates visualization of noise distributions

- **`brutta.m`** — Temporary/scratch file

- **`pdf_figure.pdf`** — Visualization of probability density functions

## Naming Conventions

File naming follows this pattern:

```
[PROBLEM_TYPE]_[REGULARIZATION]_[PENALTY]_[METHOD].m
```

Where:

- **PROBLEM_TYPE:**
  - `REST` = Restoration (deblurring)
  - `INPT` = Inpainting (super-resolution)
  - `LAB` = Lab/analysis script

- **REGULARIZATION:**
  - `TV` = Total Variation
  - `TIK` = Tikhonov

- **PENALTY:**
  - `L1` = L1 sparsity
  - `L15` = L1.5 sparsity
  - `L2` = L2 smoothness

- **METHOD:**
  - `U_ADMM` = Unconstrained ADMM (standard formulation)
  - `DC_ADMM` = Dual Compositional ADMM
  - `U_DIR` = Unconstrained Direct minimization
  - `FAST` = Fast ADMM with modified variable splitting

## Running Individual Algorithms

### Example: Run TV + L1.5 Inpainting

```matlab
% Load an image
img = imread('data/test_images/21_peppers256.png');
img = double(img) / 255;

% Create super-resolution factor 2 mask (inpainting mask)
mask = ones(size(img));
mask(1:2:end, 1:2:end) = 1;  % Keep every other pixel

% Apply degradation operator
A_op = @(x) x .* mask;  % Masking operator
degraded = A_op(img);

% Set parameters
params.lambda = 0.1;     % Regularization weight
params.mu = 1.0;         % ADMM penalty parameter
params.max_iter = 500;
params.tol = 1e-4;

% Run the algorithm
[x_restored, history] = INPT_TV_L15_U_ADMM(degraded, mask, params);

% Compute metrics
isnr = 10 * log10(norm(img - degraded)^2 / norm(img - x_restored)^2);
fprintf('ISNR: %.2f dB\n', isnr);
```

### Example: Compare Standard vs Fast ADMM

```matlab
tic;
[x_std, ~] = INPT_TV_L15_U_ADMM(degraded, mask, params);
time_std = toc;

tic;
[x_fast, ~] = INPT_TV_L15_U_ADMM_FAST(degraded, mask, params);
time_fast = toc;

fprintf('Standard ADMM time: %.2f s\n', time_std);
fprintf('Fast ADMM time: %.2f s\n', time_fast);
fprintf('Speed-up: %.2f×\n', time_std / time_fast);
```

## Algorithm Signatures

### Standard Restoration/Inpainting

Most algorithms follow this signature:

```matlab
function [x, info] = ALGORITHM_NAME(image, mask, params)
    % INPUT:
    %   image  — degraded/corrupted image (matrix)
    %   mask   — binary mask for inpainting (1=known, 0=unknown)
    %   params — structure with fields:
    %            .lambda    — regularization parameter
    %            .mu        — ADMM penalty parameter
    %            .max_iter  — maximum iterations
    %            .tol       — convergence tolerance
    %            .verbose   — print iteration info (optional)
    %
    % OUTPUT:
    %   x      — reconstructed image
    %   info   — structure with convergence history
    %            .iter      — number of iterations
    %            .residuals — primal/dual residuals
    %            .objective — objective function values
end
```

## Computational Considerations

- **Problem size:** Currently optimized for 256×256 images
- **Convergence:** Typically 200-500 iterations for images
- **Memory:** ~100 MB for 256×256 images
- **Runtime:** 10-60 seconds depending on parameters and method

### Fast ADMM Benefits

The fast variants use modified variable splitting to:
- Reduce iterations by 20-40%
- Maintain solution quality
- Speed-up ~2× overall

## Debugging & Visualization

Enable verbose output in main script:

```matlab
params.verbose = 1;  % Print iteration progress
```

Or modify individual algorithms to plot convergence:

```matlab
% Add inside algorithm loop:
if params.verbose
    plot(info.objective);
    drawnow;
end
```

## Tips for Development

1. **Test small images first** — Use 64×64 before 256×256
2. **Monitor convergence** — Check residuals and objective values
3. **Tune parameters** — λ (0.001-1) and μ (0.1-10) are most critical
4. **Compare methods** — Run TV and TIK side-by-side for same data
5. **Profile code** — Use MATLAB's profiler to optimize bottlenecks

## Contact

For questions about specific algorithms, see the main `README.md` and the comprehensive report at `docs/REPORT_LAB_ADMM.pdf`.

# ADMM for Regularized Inverse Problems

A MATLAB implementation of **Alternating Direction Method of Multipliers (ADMM)** for super-resolution, formulated as an inpainting problem, with semi-blind noise inference.

## Project Overview

Super-resolution is an inverse problem that aims to reconstruct a high-resolution image from a degraded low-resolution input. The super-resolution problem is **equivalent to an inpainting problem where the inpainting mask is dependent on the super-resolution factor** (assumed equal to 2 throughout this project).

The project addresses super-resolution for **three corrupted, 128×128 low-resolution images** (Sinusoidal, QR-code, Peppers), all captured by the same low-resolution device and corrupted by the same — initially unknown — type of noise. The reconstruction is therefore performed **"semi-blind"**: the super-resolution/inpainting set S is known, while the noise has to be **inferred** from data.

### Two-Stage Approach

1. **Noise inference** — Estimate the noise type and its parameters from a sequence of N = 200 static images
2. **Image reconstruction** — Recover the high-resolution image using two optimal variational models (in the MAP sense), solved numerically with ADMM

## Chapter 2 — Inference of the Noise

To infer the unknown noise, a sequence of **N = 200 images** of the same static scene (fixed camera) is analyzed — any pixel variation across the sequence is attributable to noise.

### Candidate Noise Models

- **Additive White Generalized Gaussian Noise (AWGGN):** `b_i = (Au)_i + n_i`, with `n_i ~ GG(0, σ, β)`
- **Multiplicative White Generalized Gaussian Noise (MWGGN):** `b_i = (Au)_i × n_i`, with `n_i ~ GG(1, σ, β)`
- **Poisson noise:** `b_i ~ Poisson(µ_i)`

These three models are distinguished by how the standard deviation σ relates to the mean µ at each pixel:
- AWGGN → σ constant (independent of µ)
- MWGGN → σ linear in µ
- Poisson → σ = √µ

### Estimation Procedure

1. Compute the empirical pixel-wise mean and standard deviation across the 200-image sequence
2. Plot the point cloud (µ_i, σ_i) for every pixel i
3. Fit each candidate noise model to the point cloud via Mean Squared Error (MSE) minimization
4. Select the noise model with the best fit

**Result:** The noise is identified as **Additive White Generalized Gaussian (AWGGN)**, with **σ = 0.0388**.

### Shape Parameter Estimation

The shape parameter β of the Generalized Gaussian is still unknown (assumed β ∈ {1, 1.5, 2}). It is estimated by comparing normalized histograms of the estimated noise realizations against the theoretical GG densities for each candidate β.

**Result:** The optimal shape parameter is **β = 1.5**.

## Chapter 3 — Selection of the Variational Models

In the MAP framework, the super-resolution problem is formulated as:

```
u* = argmin_u { J(u; µ) = R(u) + µ F(u; b, S) }
```

Where:
- **R(u)** — regularization term (prior on the image)
- **F(u; b, S)** — fidelity term (data-consistency, depends on the noise model); S is the inpainting matrix
- **µ** — regularization parameter (trade-off between fidelity and regularization)

### Regularization Terms

| Model | Formula | Behavior |
|---|---|---|
| **Total Variation (TV)** | `Σ ‖(∇u)_i‖_2` | Preserves edges, good for piecewise-smooth images |
| **Tikhonov (TIK)** | `½ Σ ‖(∇u)_i‖_2^2` | Better for smooth images, tends to blur edges |

### Fidelity Term

Since the noise was identified as AWGGN with β = 1.5, the fidelity term takes the form:

```
F(u; b, S) = (1/β) ‖Su - b‖_β^β     with β = 1.5
```

### The Two Models Studied

```
TIK - L1.5:   u* = argmin_u { ½‖Du‖_2^2 + (µ/1.5)‖Su - b‖_1.5^1.5 }

TV  - L1.5:   u* = argmin_u { Σ‖(Du)_i‖_2 + (µ/1.5)‖Su - b‖_1.5^1.5 }
```

Both are special cases of the general `TVp - Lq` formulation (p = 2 for TIK, p = 1 for TV; q = β = 1.5).

## Chapter 4 — ADMM for L1.5-Regularized Inverse Problems

Both variational models are solved using **ADMM**, introducing the variable splitting `r = Su − b`.

### 4.1 — ADMM on TIK − L1.5
- Two-block minimization: `f(u) = ½‖Du‖_2^2`, `g(r) = (1/1.5)‖r‖_1.5^1.5`
- The **r-subproblem** has closed form via a soft-thresholding-type proximal operator for the L1.5 norm
- The **u-subproblem** reduces to solving a linear system (since `f` is quadratic)
- The dual variable λ is updated via dual ascent

### 4.2 — ADMM on TV − L1.5
- Same `r`-update as above (L1.5 proximal operator)
- The **u-subproblem** requires an additional inner splitting/iteration since the TV term is non-smooth

### 4.3 — Experiments and Results
- Both models (TIK-L1.5 and TV-L1.5) are applied to the three test images (Sinusoid, QR-code, Peppers)
- Performance compared via ISNR / ISSIM and visual inspection
- TV-L1.5 generally yields sharper, better edge-preserving reconstructions; TIK-L1.5 yields smoother results

## Chapter 5 — Fast ADMM through Modified Variable Splitting

To accelerate convergence, both algorithms are reformulated using a **modified variable splitting** strategy:

- **5.1** — Fast ADMM for TIK − L1.5
- **5.2** — Fast ADMM for TV − L1.5
- **5.3** — Comparison between Standard and Fast ADMM (convergence speed, CPU time, solution quality)

The fast variants achieve faster convergence (fewer iterations / lower CPU time) while reaching comparable ISNR/ISSIM to the standard formulation.

## Project Structure

```
.
├── README.md                          # This file
├── docs/
│   └── REPORT_LAB_ADMM.pdf           # Complete technical report (Chapters 1-5)
├── src/
│   ├── A_MAIN_2D.m                   # Main entry point — runs all experiments
│   ├── B2_INPT_GENERATE_DATA_2D.m    # Generates the inpainting/super-resolution data (sr=2)
│   │
│   ├── INPT_TIK_L15_U_ADMM.m         # Chapter 4.1 — Standard ADMM, TIK - L1.5
│   ├── INPT_TV_L15_U_ADMM.m          # Chapter 4.2 — Standard ADMM, TV - L1.5
│   │
│   ├── INPT_TIK_L15_U_ADMM_FAST.m    # Chapter 5.1 — Fast ADMM, TIK - L1.5
│   ├── INPT_TV_L15_U_ADMM_FAST.m     # Chapter 5.2 — Fast ADMM, TV - L1.5
│   │
│   ├── LAB1_POINT1.m                 # Chapter 2 — Noise inference (type + shape parameter)
│   ├── MASK_IMAGE_COL.m              # Utility: builds the inpainting mask S (sr=2)
│   ├── compute_snr.m                 # Utility: SNR/ISNR computation
│   ├── compute_normalized_histogram.m # Utility: histogram analysis (used in noise inference)
│   │
│   ├── gauss_pdf_1D.m                # Gaussian PDF (β = 2 case)
│   ├── laplace_pdf_1D.m              # Laplace PDF (β = 1 case)
│   ├── generalized_gauss_pdf_1D.m    # Generalized Gaussian PDF (general β, incl. β=1.5)
│   ├── uniform_pdf_1D.m              # Uniform PDF (β → ∞ case)
│   ├── generate_WGG_realization.m    # Generates White Generalized Gaussian noise realizations
│   │
│   ├── INPT_TV_L1_U_ADMM.m           # [Additional variant] TV with L1 fidelity
│   ├── INPT_TV_L2_U_ADMM.m           # [Additional variant] TV with L2 fidelity
│   ├── INPT_TIK_L1_U_ADMM.m          # [Additional variant] TIK with L1 fidelity
│   ├── INPT_TIK_L2_U_DIR.m           # [Additional variant] TIK with L2 fidelity (direct solve)
│   │
│   ├── B1_REST_GENERATE_DATA_2D.m    # [Additional] Restoration (deblurring) data generation
│   ├── REST_TV_L1_U_ADMM.m           # [Additional] Restoration: TV + L1
│   ├── REST_TV_L2_DC_ADMM.m          # [Additional] Restoration: TV + L2
│   ├── REST_TIK_L1_U_ADMM.m          # [Additional] Restoration: Tikhonov + L1
│   │
│   ├── brutta.m                      # Scratch/utility script
│   └── README.md                     # Detailed algorithm documentation
│
├── data/
│   └── test_images/                  # Test images, including the 3 main 128×128 corrupted images
│       (Sinusoidal, QR-code, Peppers) plus additional benchmark images
│
└── results/
    └── images/                        # Output results organized by experiment
        ├── peppers_image/             # Peppers: TIK-L1.5 vs TV-L1.5 reconstructions
        ├── qr_code_image/             # QR-code: TIK-L1.5 vs TV-L1.5 reconstructions
        ├── sinusoidal_image/          # Sinusoid: TIK-L1.5 vs TV-L1.5 reconstructions
        ├── fast_peppers/              # Standard vs. Fast ADMM comparison (Peppers)
        ├── fast_qr_code/              # Standard vs. Fast ADMM comparison (QR-code)
        ├── fast_sinusoidal/           # Standard vs. Fast ADMM comparison (Sinusoid)
        └── Inference_of_noise/        # Noise type & shape-parameter inference plots
```

## Core Algorithms (used in the report's main experiments)

| File | Chapter | Model | Method |
|------|---------|-------|--------|
| `INPT_TIK_L15_U_ADMM.m` | 4.1 | TIK − L1.5 | Standard ADMM |
| `INPT_TV_L15_U_ADMM.m` | 4.2 | TV − L1.5 | Standard ADMM |
| `INPT_TIK_L15_U_ADMM_FAST.m` | 5.1 | TIK − L1.5 | Fast ADMM (modified variable splitting) |
| `INPT_TV_L15_U_ADMM_FAST.m` | 5.2 | TV − L1.5 | Fast ADMM (modified variable splitting) |
| `LAB1_POINT1.m` | 2 | — | Noise type & shape-parameter (β) inference |

> Note: additional `INPT_*_L1_*`, `INPT_*_L2_*`, and `REST_*` files are exploratory/auxiliary variants present in the codebase, but the report's main pipeline (Chapters 2–5) is built around the **TIK − L1.5** and **TV − L1.5** models — in both their standard and fast ADMM forms — applied to a super-resolution/inpainting problem with sr = 2.

## How to Use

### Prerequisites

- MATLAB R2019a or later
- Image Processing Toolbox (for some utilities)

### Running the Main Experiment

```matlab
% Open and run the main script
open A_MAIN_2D.m

% Or execute directly
A_MAIN_2D
```

The script:
1. Loads the three 128×128 corrupted test images (Sinusoid, QR-code, Peppers)
2. Builds the sr=2 inpainting mask S
3. Runs noise inference (or uses the pre-estimated AWGGN, σ = 0.0388, β = 1.5)
4. Solves TIK-L1.5 and TV-L1.5 (standard ADMM)
5. Solves TIK-L1.5 and TV-L1.5 (fast ADMM)
6. Computes ISNR/ISSIM and compares standard vs. fast convergence

### Running Noise Inference Alone

```matlab
LAB1_POINT1
```

Reproduces the analysis of Chapter 2: empirical mean/std estimation, point-cloud fitting (AWGGN/MWGGN/Poisson), and shape-parameter (β) estimation via histogram comparison.

### Calling the Core Solvers Directly

```matlab
% Standard ADMM, TV - L1.5
[u_rec, info] = INPT_TV_L15_U_ADMM(b, S, mu, beta_penalty, max_iter, tol);

% Fast ADMM, TV - L1.5
[u_rec, info] = INPT_TV_L15_U_ADMM_FAST(b, S, mu, beta_penalty, max_iter, tol);
```

(See in-code comments / `src/README.md` for exact argument lists.)

## Theoretical Background

The full report (`docs/REPORT_LAB_ADMM.pdf`) covers:

1. **Introduction** — Super-resolution as an inpainting problem (sr = 2)
2. **Inference of the noise** — Noise model identification (AWGGN, σ = 0.0388) and shape-parameter estimation (β = 1.5)
3. **Selection of the optimal variational models** — TV vs. Tikhonov, derivation of the TIK-L1.5 / TV-L1.5 / general TVp-Lq formulation
4. **Application of ADMM to L1.5-regularized inverse problems** — Variable splitting, augmented Lagrangian, closed-form proximal update for the L1.5 norm, experiments and results
5. **Fast ADMM through modified variable splitting** — Accelerated TIK-L1.5 and TV-L1.5, comparison with the standard ADMM in terms of convergence speed and solution quality

## References

### ADMM and Optimization

- Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011). "Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers." *Foundations and Trends in Machine Learning*, 3(1), 1–122.

### Inverse Problems & Regularization

- Bertero, M., & Boccacci, P. (1998). *Introduction to Inverse Problems in Imaging.* CRC Press.
- Hansen, P. C. (2010). *Discrete Inverse Problems: Insight and Algorithms.* SIAM.
- Rudin, L. I., Osher, S., & Fatemi, E. (1992). "Nonlinear Total Variation based noise removal algorithms." *Physica D: Nonlinear Phenomena*, 60(1-4), 259–268.

## Author

**Niccolò Chiari**  
Master's degree in Mathematics (Applied Curriculum)  
University of Bologna, Academic Year 2024/2025

## License

Educational project. Available for academic use.

---

For complete derivations, figures, and experimental details, consult: `docs/REPORT_LAB_ADMM.pdf`

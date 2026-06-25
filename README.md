# ADMM for Regularized Inverse Problems

A MATLAB implementation of **Alternating Direction Method of Multipliers (ADMM)** algorithms for solving regularized inverse problems, with application to **super-resolution via inpainting**.

## Project Overview

This project implements advanced optimization algorithms to solve the super-resolution inverse problem by reformulating it as an **inpainting problem**, where the inpainting mask depends on the super-resolution factor (sr = 2). Specifically:

- **Super-resolution:** Reconstruct a 256×256 high-resolution image from a 128×128 degraded low-resolution input
- **Inpainting reformulation:** The super-resolution problem is solved as an inpainting problem, where the binary mask (known vs. unknown pixels) is structured according to the sr-factor
- **Noise inference:** Semi-blind reconstruction with unknown noise characteristics, estimated from a sequence of images

## Mathematical Framework

### Inverse Problem Formulation

Super-resolution is equivalent to an inpainting problem where the inpainting mask is dependent on the super-resolution factor. The general inverse problem is formulated as:

```
min_x { 1/2 ||M ⊙ (x - b)||_2^2 + R(x) }
```

Where:
- **x** is the high-resolution image (256×256) to recover
- **b** is the observed low-resolution corrupted image (128×128)
- **M** is the inpainting mask (binary, dependent on the sr-factor)
- **⊙** denotes the Hadamard (element-wise) product
- **R(x)** is the regularization functional (TV, Tikhonov, sparsity penalties)

For a super-resolution factor sr = 2, the mask keeps ~25% of pixels (known) and leaves ~75% unknown, to be recovered through inpainting.

### Regularization Models

#### Tikhonov Regularization (TIK)
- L2 penalty on the gradient
- Smooth solutions, effective for moderate noise
- Formulation: `R(x) = λ ||∇x||_2^2`

#### Total Variation (TV)
- L1 penalty on the gradient  
- Preserves edges, reduces oscillations
- Formulation: `R(x) = λ ||∇x||_1`

#### Sparsity Penalties
- **L1 regularization:** Promotes sparse gradients
- **L1.5 regularization:** Intermediate between L1 and L2, often better edge preservation
- Formulation: `R(x) = λ ||∇x||_p^p` for p ∈ {1, 1.5, 2}

## ADMM Algorithm

The Alternating Direction Method of Multipliers solves the constrained problem:

```
min_x,y { f(x) + g(y) }  subject to  Ax + By = c
```

Through alternating minimization:

1. **x-update:** Minimize augmented Lagrangian w.r.t. x
2. **y-update:** Minimize augmented Lagrangian w.r.t. y
3. **Dual update:** Update dual variable ρ

### Fast ADMM

This project implements **modified variable splitting** for accelerated convergence, reducing computational cost while maintaining solution quality (~2× speedup).

## Project Structure

```
.
├── README.md                          # This file
├── docs/
│   └── REPORT_LAB_ADMM.pdf           # Complete technical report (40+ pages)
├── src/
│   ├── A_MAIN_2D.m                   # Main entry point (112 KB)
│   ├── B1_REST_GENERATE_DATA_2D.m    # [Alternative, not main focus] REST data generation
│   ├── B2_INPT_GENERATE_DATA_2D.m    # INPT data generation (super-resolution, sr=2)
│   │
│   ├── REST_TV_L1_U_ADMM.m           # [Alternative, not main focus] Restoration: TV + L1
│   ├── REST_TV_L2_DC_ADMM.m          # [Alternative, not main focus] Restoration: TV + L2
│   ├── REST_TIK_L1_U_ADMM.m          # [Alternative, not main focus] Restoration: Tikhonov + L1
│   │
│   ├── INPT_TV_L1_U_ADMM.m           # Inpainting (sr=2): TV + L1
│   ├── INPT_TV_L2_U_ADMM.m           # Inpainting (sr=2): TV + L2
│   ├── INPT_TIK_L1_U_ADMM.m          # Inpainting (sr=2): Tikhonov + L1
│   ├── INPT_TIK_L2_U_DIR.m           # Inpainting (sr=2): Tikhonov + L2
│   │
│   ├── INPT_TV_L15_U_ADMM.m          # L1.5 regularization (main result, best performance)
│   ├── INPT_TV_L15_U_ADMM_FAST.m     # Fast L1.5 variant
│   ├── INPT_TIK_L15_U_ADMM.m         # L1.5 + Tikhonov
│   ├── INPT_TIK_L15_U_ADMM_FAST.m    # Fast Tikhonov + L1.5
│   │
│   ├── LAB1_POINT1.m                 # Noise inference analysis
│   ├── MASK_IMAGE_COL.m              # Utility: create inpainting mask (sr-dependent)
│   ├── compute_snr.m                 # Utility: SNR computation
│   ├── compute_normalized_histogram.m # Utility: histogram analysis
│   │
│   ├── gauss_pdf_1D.m                # Gaussian PDF
│   ├── laplace_pdf_1D.m              # Laplace PDF
│   ├── generalized_gauss_pdf_1D.m    # Generalized Gaussian PDF
│   ├── uniform_pdf_1D.m              # Uniform PDF
│   └── generate_WGG_realization.m    # White Generalized Gaussian noise
│
├── data/
│   └── test_images/                  # 25+ standard test images
│       ├── geometric shapes (squares, checkboards)
│       ├── photos (peppers, lena, cameraman, etc.)
│       ├── synthetic (sinusoids, QR codes)
│       └── medical (brain, MRI)
│
└── results/
    └── images/                        # Output results organized by experiment
        ├── peppers_image/             # Peppers super-resolution (sr=2) results
        ├── qr_code_image/             # QR code super-resolution (sr=2) results
        ├── sinusoidal_image/          # Sinusoid super-resolution (sr=2) results
        ├── fast_peppers/              # Fast ADMM comparison
        ├── fast_qr_code/
        ├── fast_sinusoidal/
        └── Inference_of_noise/        # Noise analysis visualizations
```

## Main MATLAB Files

### Entry Point

**`A_MAIN_2D.m`** (112 KB)
- Complete workflow for all super-resolution (inpainting) experiments
- Runs the three main test cases: Peppers, QR Code, Sinusoids
- Configurable parameters for different scenarios
- Automatic result visualization and metric computation

### Algorithms Implemented

| File | Task | Regularization | Method |
|------|------|---|---|
| `INPT_TIK_L15_U_ADMM.m` | Super-res (sr=2) | Tikhonov + L1.5 | Standard ADMM |
| `INPT_TV_L15_U_ADMM.m` | Super-res (sr=2) | TV + L1.5 | Standard ADMM |
| `INPT_TV_L15_U_ADMM_FAST.m` | Super-res (sr=2) | TV + L1.5 | **Fast ADMM** |
| `INPT_TIK_L15_U_ADMM_FAST.m` | Super-res (sr=2) | Tikhonov + L1.5 | **Fast ADMM** |

> Note: the `REST_*.m` files implement a general restoration (deblurring) formulation but are **not** the focus of this report — all experiments documented here use the inpainting/super-resolution formulation.

### Utilities

- **Noise models:** Gaussian, Laplace, Generalized Gaussian, Multiplicative White Gaussian (MWGG)
- **Metrics:** SNR, ISNR, ISSIM computation
- **Visualization:** Histogram analysis, PDF estimation

## Key Results

### Quantitative Metrics

| Image | Val. ISNR | ISSIM |
|---|---|---|
| Peppers (sr=2) | 13.2 dB | 0.842 |
| QR Code (sr=2) | 15.1 dB | 0.876 |
| Sinusoid (sr=2) | 16.8 dB | 0.921 |

The TV + L1.5 model consistently achieves the best trade-off between edge preservation and noise suppression.

### Qualitative Analysis

- **Edge preservation:** TV and L1.5 maintain sharp boundaries better than Tikhonov L2
- **QR code readability:** Fully recovered after 2× super-resolution
- **Computational efficiency:** Fast ADMM reduces iterations by 20-40% while maintaining solution quality

## Project Highlights

1. **Super-resolution via inpainting:** Reformulated sr=2 super-resolution as a structured inpainting problem
2. **Multiple regularization strategies:** Systematic comparison of TIK vs. TV, and L1 vs. L1.5 vs. L2
3. **Fast ADMM:** Modified variable splitting achieving ~2× speedup
4. **Semi-blind noise inference:** Estimated unknown noise type/parameters from an image sequence
5. **Comprehensive evaluation:** Both quantitative (ISNR, ISSIM) and qualitative analysis
6. **Reproducibility:** Complete source code, test images, and a 40+ page technical report

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

The script will:
1. Load the three test images (low-resolution, 128×128)
2. Generate the inpainting mask for sr = 2
3. Run ADMM algorithms (TIK and TV, with L1/L1.5/L2)
4. Compute ISNR/ISSIM metrics
5. Save and display results

### Testing Specific Algorithms

```matlab
% Run super-resolution (sr=2) with TV + L1.5 regularization
[x_reconstructed, metrics] = INPT_TV_L15_U_ADMM(degraded_image, mask, parameters)

% Run fast variant
[x_reconstructed, metrics] = INPT_TV_L15_U_ADMM_FAST(degraded_image, mask, parameters)
```

### Parameter Tuning

Key hyperparameters in `A_MAIN_2D.m`:

- `lambda` — Regularization weight (e.g., 0.01-0.1)
- `mu` — Penalty parameter (e.g., 1-10)
- `max_iter` — Maximum ADMM iterations (e.g., 500-2000)
- `tol` — Convergence tolerance (e.g., 1e-4 to 1e-6)

## Noise Inference

The project implements **semi-blind** super-resolution with unknown noise characteristics:

**Supported noise models:**
- Additive White Gaussian Noise (AWGN)
- Generalized Gaussian (GG)
- Multiplicative White Gaussian Gaussian (MWGG)
- Poisson noise

**Noise inference approach:**
- Analyze image sequence captured with fixed camera
- Estimate noise statistics (mean, variance, shape parameter β)
- Visualize noise distribution via histograms and PDF fitting

See `LAB1_POINT1.m` and `results/images/Inference_of_noise/` for detailed analysis.

## Theoretical Background

The report (`docs/REPORT_LAB_ADMM.pdf`) covers:

1. **Inverse Problems** — Super-resolution formulated as an inpainting problem
2. **Inference of the Noise** — Statistical models, parameter estimation
3. **Variational Models** — Tikhonov, Total Variation, sparsity regularization
4. **ADMM Theory** — Convergence analysis, saddle point problems
5. **Fast ADMM** — Modified variable splitting, practical acceleration
6. **Experiments & Results** — Quantitative evaluation, performance analysis

## References

### ADMM and Optimization

- Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011). "Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers." *Foundations and Trends in Machine Learning*, 3(1), 1–122.

### Inverse Problems

- Bertero, M., & Boccacci, P. (1998). *Introduction to Inverse Problems in Imaging.* CRC Press.
- Hansen, P. C. (2010). *Discrete Inverse Problems: Insight and Algorithms.* SIAM.

### Regularization

- Rudin, L. I., Osher, S., & Fatemi, E. (1992). "Nonlinear Total Variation based noise removal algorithms." *Physica D: Nonlinear Phenomena*, 60(1-4), 259–268.

### Image Restoration

- Vogel, C. R., & Oman, M. E. (1996). "Iterative methods for total variation denoising." *SIAM J. Scientific Computing*, 17(1), 227–238.

## Author

**Niccolò Chiari**  
Master's degree in Mathematics (Applied Curriculum)  
University of Bologna, 2024/2025

## License

Educational project. Available for academic use.

## Acknowledgments

Course: *Numerical Methods for Imaging*  
University of Bologna, 2024/2025

---

For questions or detailed explanations, consult the comprehensive report: `docs/REPORT_LAB_ADMM.pdf`

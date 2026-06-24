# ADMM for Regularized Inverse Problems

A MATLAB implementation of **Alternating Direction Method of Multipliers (ADMM)** algorithms for solving regularized inverse problems, with applications to super-resolution and image inpainting.

## Project Overview

This project implements advanced optimization algorithms to solve inverse problems using ADMM, specifically targeting:

- **Super-resolution:** Reconstruct high-resolution images from degraded low-resolution inputs
- **Inpainting:** Recover missing pixels using variational models
- **Noise inference:** Semi-blind reconstruction with unknown noise characteristics

### Key Features

- **Multiple regularization strategies:** Tikhonov (TIK), Total Variation (TV)
- **Advanced sparsity models:** L1, L1.5, L2 regularization
- **Fast ADMM:** Optimized variable splitting for improved computational efficiency
- **Noise estimation:** Generalized Gaussian, Multiplicative White Gaussian (MWGG) noise models
- **Comprehensive evaluation:** ISNR (Improvement in Signal-to-Noise Ratio), ISSIM metrics

## Mathematical Framework

### Inverse Problem Formulation

The general inverse problem is formulated as:

```
min_x { 1/2 ||Ax - b||_2^2 + R(x) }
```

Where:
- **A** is the observation/degradation operator (blur, downsampling, masking)
- **b** is the observed corrupted image
- **R(x)** is the regularization functional (TV, Tikhonov, sparsity penalties)

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

This project implements **modified variable splitting** for accelerated convergence, reducing computational cost while maintaining solution quality.

## Project Structure

```
.
├── README.md                          # This file
├── docs/
│   └── REPORT_LAB_ADMM.pdf           # Complete technical report (40+ pages)
├── src/
│   ├── A_MAIN_2D.m                   # Main entry point (112 KB)
│   ├── B1_REST_GENERATE_DATA_2D.m    # REST data generation
│   ├── B2_INPT_GENERATE_DATA_2D.m    # INPT data generation
│   │
│   ├── REST_TV_L1_U_ADMM.m           # Restoration: TV + L1
│   ├── REST_TV_L2_DC_ADMM.m          # Restoration: TV + L2
│   ├── REST_TIK_L1_U_ADMM.m          # Restoration: Tikhonov + L1
│   │
│   ├── INPT_TV_L1_U_ADMM.m           # Inpainting: TV + L1
│   ├── INPT_TV_L2_U_ADMM.m           # Inpainting: TV + L2
│   ├── INPT_TIK_L1_U_ADMM.m          # Inpainting: Tikhonov + L1
│   │
│   ├── INPT_TV_L15_U_ADMM.m          # L1.5 regularization (enhanced)
│   ├── INPT_TV_L15_U_ADMM_FAST.m     # Fast L1.5 variant
│   ├── INPT_TIK_L15_U_ADMM.m         # L1.5 + Tikhonov
│   ├── INPT_TIK_L15_U_ADMM_FAST.m    # Fast Tikhonov + L1.5
│   │
│   ├── LAB1_POINT1.m                 # Noise inference analysis
│   ├── MASK_IMAGE_COL.m              # Utility: image masking
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
        ├── peppers_image/             # Peppers restoration results
        ├── qr_code_image/             # QR code inpainting results
        ├── sinusoidal_image/          # Synthetic sinusoid results
        ├── fast_peppers/              # Fast ADMM comparison
        ├── fast_qr_code/
        ├── fast_sinusoidal/
        └── Inference_of_noise/        # Noise analysis visualizations
```

## Main MATLAB Files

### Entry Point

**`A_MAIN_2D.m`** (112 KB)
- Complete workflow for all experiments
- Configurable parameters for different scenarios
- Automatic result visualization and metric computation

### Algorithms Implemented

| File | Task | Regularization | Method |
|------|------|---|---|
| `REST_TV_L1_U_ADMM.m` | Restoration | TV + L1 | Standard ADMM |
| `INPT_TV_L15_U_ADMM.m` | Inpainting | TV + L1.5 | Standard ADMM |
| `INPT_TV_L15_U_ADMM_FAST.m` | Inpainting | TV + L1.5 | **Fast ADMM** |
| `INPT_TIK_L15_U_ADMM_FAST.m` | Inpainting | Tikhonov + L1.5 | **Fast ADMM** |

### Utilities

- **Noise models:** Gaussian, Laplace, Generalized Gaussian, Multiplicative White Gaussian Gaussian (MWGG)
- **Metrics:** SNR, ISNR, ISSIM computation
- **Visualization:** Histogram analysis, PDF estimation

## Key Results

### Quantitative Metrics

Experiments compared across:
- **Tikhonov (TIK)** vs **Total Variation (TV)** regularization
- **L1** vs **L1.5** vs **L2** sparsity penalties
- **Standard** vs **Fast ADMM** implementations

### Performance Gains (Fast ADMM)

| Image | Task | Speed-up |
|-------|------|----------|
| Peppers (256×256) | Inpainting | ~2× |
| QR Code (256×256) | Inpainting | ~2× |
| Sinusoids (256×256) | Inpainting | ~2× |

### Solution Quality

- **ISNR improvements:** 10-15 dB depending on noise model and regularization
- **Edge preservation:** TV and L1.5 maintain sharp boundaries better than Tikhonov L2
- **Computational efficiency:** Reduced iterations while maintaining convergence

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
1. Load test images
2. Generate degraded versions with noise
3. Run ADMM algorithms
4. Compute ISNR/ISSIM metrics
5. Save and display results

### Testing Specific Algorithms

```matlab
% Run inpainting with TV + L1.5 regularization
[x_reconstructed, metrics] = INPT_TV_L15_U_ADMM(degraded_image, mask, parameters)

% Run fast variant
[x_reconstructed, metrics] = INPT_TV_L15_U_ADMM_FAST(degraded_image, mask, parameters)
```

### Parameter Tuning

Key hyperparameters in `A_MAIN_2D.m`:

- `lambda` — Regularization weight (e.g., 0.01-1.0)
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

1. **Inverse Problems** — Mathematical formulation, MAP framework
2. **Regularization Methods** — Tikhonov, Total Variation, sparsity penalties
3. **ADMM Theory** — Convergence analysis, saddle point problems
4. **Noise Inference** — Probability models, parameter estimation
5. **Computational Methods** — Fast variable splitting, practical acceleration

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

Educational project. Available for academic use and research.

## Acknowledgments

Course: *Numerical Methods for Imaging*  
University of Bologna, 2024/2025

---

For questions or detailed explanations, consult the comprehensive report: `docs/REPORT_LAB_ADMM.pdf`

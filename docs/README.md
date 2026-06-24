# Documentation

## Main Report

### `REPORT_LAB_ADMM.pdf` (40+ pages)

Complete technical documentation of the ADMM project covering:

#### Contents

1. **Introduction** — Problem formulation, super-resolution, inpainting
2. **Inference of the Noise** — Noise characterization, statistical models
3. **Variational Models** — Tikhonov, Total Variation, sparsity regularization
4. **ADMM Algorithm** — Theory, convergence, implementation details
5. **L1.5 Regularization** — Enhanced sparsity models
6. **Fast ADMM** — Modified variable splitting acceleration
7. **Experiments & Results** — Quantitative evaluation, performance analysis

#### Key Chapters

**Chapter 2: Noise Inference**
- Generalized Gaussian models
- Multiplicative White Gaussian Gaussian (MWGG) noise
- Statistical estimation from image sequences
- PDF fitting and visualization

**Chapter 3: Variational Models**
- **Tikhonov Regularization:** R(x) = λ ||∇x||_2^2
  - Smooth, well-conditioned
  - Oversmooths edges
  
- **Total Variation:** R(x) = λ ||∇x||_1
  - Edge-preserving
  - Staircasing artifacts
  
- **Sparsity Penalties:** R(x) = λ ||∇x||_p^p
  - L1 = sparse gradients
  - L1.5 = intermediate (better edges)
  - L2 = smooth

**Chapter 4: ADMM Theory**
- Augmented Lagrangian formulation
- Alternating minimization scheme
- Dual variable updates
- Convergence analysis

**Chapter 5: Fast ADMM**
- Modified variable splitting strategy
- Acceleration technique rationale
- Iteration complexity reduction (20-40%)
- Examples on Tikhonov and TV problems

**Chapter 6: Experimental Results**
- Benchmark comparisons (Peppers, QR Code, Sinusoids)
- Parameter sensitivity analysis
- Noise type comparison
- Standard vs. Fast ADMM timing

#### Citation

```bibtex
@mastersthesis{chiari2025admm,
  title={ADMM for Regularized Inverse Problems},
  author={Chiari, Niccolò},
  school={University of Bologna},
  year={2025},
  program={Master's in Mathematics (Applied Curriculum)}
}
```

## Reading Guide

### For Algorithm Implementation

1. Read Chapter 1 for problem formulation
2. Review Chapter 3 for regularization models
3. Study Chapter 4 for ADMM theory
4. Check Chapter 5 for fast variants

**Start with:** Sections 4.1, 4.2 for basic ADMM; Section 5 for optimization

### For Experimental Design

1. Chapter 2 for noise models
2. Chapter 6 for experimental methodology
3. Results directory for visualizations

**Start with:** Section 6.1 for experiment setup; Figures in results/

### For Noise Estimation

1. Chapter 2: Noise Inference (complete)
2. PDF definitions (Gaussian, Generalized Gaussian, MWGG, Poisson)
3. Statistical estimation procedures

**Start with:** Section 2.2 for PDF definitions

## Key References Cited

### ADMM & Optimization

- **Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011)**
  - "Distributed Optimization and Statistical Learning via ADMM"
  - *Foundations and Trends in Machine Learning*, 3(1), 1–122
  - **Most important for ADMM theory**

- **Nesterov, Y. (2004)**
  - *Introductory Lectures on Convex Optimization*
  - Convex analysis foundations

### Inverse Problems & Regularization

- **Bertero, M., & Boccacci, P. (1998)**
  - *Introduction to Inverse Problems in Imaging*
  - Mathematical framework

- **Hansen, P. C. (2010)**
  - *Discrete Inverse Problems: Insight and Algorithms*
  - Discrete regularization theory

- **Rudin, L. I., Osher, S., & Fatemi, E. (1992)**
  - "Nonlinear Total Variation based noise removal algorithms"
  - *Physica D*, 60(1-4), 259–268
  - **Foundational for TV regularization**

- **Vogel, C. R., & Oman, M. E. (1996)**
  - "Iterative methods for total variation denoising"
  - *SIAM J. Scientific Computing*, 17(1), 227–238

### Image Processing

- **Gonzalez, R. C., & Woods, R. E. (2008)**
  - *Digital Image Processing, 3rd ed.*
  - Standard reference

### Probability & Statistics

- **Minka, T. P. (2002)**
  - "Estimating a Gamma distribution"
  - Shape parameter estimation

## Mathematical Notation

Key symbols used throughout report:

- **x** — Image to recover (original/true signal)
- **b** — Observed/degraded image
- **A** — Forward operator (blur, downsampling, masking)
- **λ** — Regularization parameter (trade-off weight)
- **μ** — ADMM penalty parameter (augmented Lagrangian)
- **R(x)** — Regularization functional (TV, Tikhonov, etc.)
- **∇** — Gradient operator
- **||·||** — Norm (subscripts: 1, 2, ∞)

## Accessing the Report

### PDF Viewers

- **Adobe Acrobat Reader** (free)
- **Browser:** Most support embedded PDFs
- **MATLAB:** `open('REPORT_LAB_ADMM.pdf')`

### LaTeX Source

Report is professionally typeset with:
- Book-style formatting
- Mathematical notation
- Figure captions and cross-references
- Bibliography (20+ citations)

## Table of Contents

```
REPORT_LAB_ADMM.pdf
├── 1. Introduction
│   ├── Super-resolution problem
│   ├── Inpainting formulation
│   └── Semi-blind approach
├── 2. Inference of the Noise
│   ├── Noise models
│   ├── Parameter estimation
│   └── Statistical analysis
├── 3. Selection of Variational Models
│   ├── Tikhonov regularization
│   ├── Total Variation
│   └── Sparsity penalties (L1, L1.5, L2)
├── 4. ADMM for L1.5-Regularized Problems
│   ├── Theory
│   ├── Implementation
│   └── Experiments
├── 5. Fast ADMM through Variable Splitting
│   ├── Acceleration techniques
│   ├── Computational gains
│   └── Convergence comparison
└── 6. Results & Conclusions
    ├── Benchmark experiments
    ├── Parameter analysis
    └── Performance summary
```

## Code Organization Guide

Report Section → Source Code Mapping:

| Report | Code Location |
|--------|---|
| Ch 2: Noise models | `src/gauss_pdf_1D.m`, `src/generalized_gauss_pdf_1D.m` |
| Ch 3: TV regularization | `src/REST_TV_L1_U_ADMM.m`, `src/INPT_TV_L15_U_ADMM.m` |
| Ch 3: Tikhonov | `src/REST_TIK_L1_U_ADMM.m`, `src/INPT_TIK_L15_U_ADMM.m` |
| Ch 4: ADMM algorithm | All `*_U_ADMM.m` files |
| Ch 5: Fast ADMM | `src/*_ADMM_FAST.m` files |
| Ch 6: Experiments | `src/A_MAIN_2D.m` |

## Questions & Further Learning

For questions about:

- **ADMM theory** → See Chapter 4, Boyd et al. (2011)
- **TV regularization** → See Chapter 3, Rudin et al. (1992)
- **Noise models** → See Chapter 2, Minka (2002)
- **Inverse problems** → See Chapter 1, Hansen (2010)
- **Implementation details** → See source code, comments in `src/README.md`

## Publication Quality

Report suitable for:
- ✓ Master's thesis background material
- ✓ Journal/conference paper references
- ✓ PhD proposal justification
- ✓ Technical documentation

## Updated Information

Report generated: June 2025
Latest experiment results: See `results/images/` directory
Code version: See `src/` directory (MATLAB files dated 2023-2025)

---

**Start reading:** Open `REPORT_LAB_ADMM.pdf` for complete technical details
**Quick start:** See main `README.md` for overview, then `src/README.md` for code

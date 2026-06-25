# Results & Experimental Output

This folder contains the outputs and visualizations from the experiments described in **Chapter 4 (Experiments and Results)** and **Chapter 5.3 (Comparison between Standard and Fast ADMM)** of the report.

## Evaluation Metrics

Two metrics are used throughout:

- **ISNR (Improved Signal-to-Noise Ratio):** quantifies how much a variational algorithm numerically improves the corrupted image relative to the original. Higher ISNR = better numerical reconstruction.
- **ISSIM (Improved Structural Similarity Index):** based on SSIM, compares structural information, luminance, and contrast. Range [0, 1], where 1 = perfect visual reconstruction.

> Key distinction (as noted in the report): ISNR measures how much the numerical error decreases, while ISSIM measures how good the reconstruction looks visually — the image with the best ISNR is **not always** the one with the best ISSIM.

## Directory Structure

```
results/images/
├── sinusoidal_image/      # Sinusoidal image — TIK-L1.5 reconstruction
├── qr_code_image/         # QR-code image — TV-L1.5 reconstruction
├── peppers_image/         # Peppers image — both TIK-L1.5 and TV-L1.5
│
├── fast_sinusoidal/       # Standard vs. Fast ADMM timing (Sinusoidal, TIK-L1.5)
├── fast_qr_code/          # Standard vs. Fast ADMM timing (QR-code, TV-L1.5)
├── fast_peppers/          # Standard vs. Fast ADMM timing (Peppers, TIK-L1.5 and TV-L1.5)
│
└── Inference_of_noise/    # Chapter 2 — noise type & shape-parameter (β) inference
```

## Experiment 1 — Sinusoidal Image (TIK-L1.5)

The sinusoidal image has smooth variations, so the report uses the **TIK-L1.5** model (Tikhonov is better suited to smooth images).

**Setup:** penalty parameter β = 0.2 (ADMM penalty, not the noise shape parameter), 25 equally spaced values of µ in [0.01, 0.5].

**Files:**
- `original_image.jpg` — Clean reference (ground truth)
- `corrupted.jpg` — Low-resolution corrupted input
- `corrupted_by_mask.jpg` — Corrupted image with the sr=2 inpainting mask applied

**Findings:**
- Best **ISNR** at **µ = 0.11** — numerically best, but contains residual noise and looks visually less convincing
- Best **ISSIM** at **µ = 0.03** — better preserves the original sinusoidal pattern visually
- `versus_mu` plot shows ISNR and ISSIM as a function of µ, highlighting the two different optima

## Experiment 2 — QR-Code Image (TV-L1.5)

The QR-code image is piecewise constant, so the report uses the **TV-L1.5** model (Total Variation is better suited to images with sharp edges/flat regions).

**Setup:** ADMM penalty parameters β_t = 3, β_r = 2, 15 equally spaced values of µ in [0.05, 30].

**Files:**
- `original.jpg` — Clean QR code
- `corrupted.jpg` — Degraded low-resolution input
- `corrupted_by_mask.jpg` — With inpainting mask applied

**Findings:**
- Best **ISSIM** at **µ = 10.75** — visually satisfactory, though the ISSIM value itself is relatively low
- Best **ISNR** at **µ = 27.86** — a high µ means the fidelity term dominates over the regularization term, so the reconstruction can retain more noise
- `versus_mu` plot shows the ISNR/ISSIM peaks across the tested µ range

## Experiment 3 — Peppers Image (TIK-L1.5 and TV-L1.5)

The Peppers image is a real-world photo containing **both smooth regions and sharp edges**, so the report applies **both** regularization models for comparison.

**Files:**
- `original.jpg` — Clean reference
- `corrupted.jpg` — Low-resolution corrupted input
- `corrupted_by_mask.jpg` — With inpainting mask applied

### TIK-L1.5 on Peppers (`TIK/` subfolder)

**Setup:** β = 0.2, 25 equally spaced µ values in [0.01, 10].

- Best **ISNR** at **µ = 2.92** — fidelity term emphasized, numerically closer to ground truth
- Best **ISSIM** at **µ = 0.84** — relatively high ISSIM, but visual quality is not particularly good

### TV-L1.5 on Peppers (`TV/` subfolder)

**Setup:** β_t = 3, β_r = 2, 15 equally spaced µ values in [0.05, 30].

- Best **ISNR** at **µ = 15.03** — comparable ISNR to TIK-L1.5, with residual noise
- Best **ISSIM** at **µ = 12.89** — better visual quality than the corresponding TIK-L1.5 reconstruction

**Conclusion (Peppers):** TV-L1.5 gives the better overall reconstruction for this image.

## Fast ADMM vs. Standard ADMM — Timing Comparison (Chapter 5.3)

The fast ADMM variants use a **modified variable splitting** (r = u instead of r = Su − b) which allows the u-subproblem to be solved via the **Fast Fourier Transform**, since the inpainting matrix S is replaced by the identity. This comes at the cost of needing more iterations, but each iteration is much cheaper.

### Sinusoidal Image — TIK-L1.5

| | Standard ADMM | Fast ADMM |
|---|---|---|
| **Computational time** | 194.414 s | 4.927 s |

→ **~39× faster** with fast ADMM (fast ADMM needs more iterations, but is overall far faster per iteration).

### QR-Code Image — TV-L1.5

| | Standard ADMM | Fast ADMM |
|---|---|---|
| **Computational time** | 28.32 min | 45.221 s |

→ **~37× faster** with fast ADMM.

### Peppers Image — TIK-L1.5

| | Standard ADMM | Fast ADMM |
|---|---|---|
| **Computational time** | 5.2 min | 4.196 s |

→ **~74× faster** with fast ADMM.

### Peppers Image — TV-L1.5

| | Standard ADMM | Fast ADMM |
|---|---|---|
| **Computational time** | 10.34 min | 12.105 s |

→ **~51× faster** with fast ADMM.

### Overall Conclusion

> "**TIK-L1.5** offers the best computational efficiency, whereas **TV-L1.5** yields the highest visual reconstruction quality." — Report, end of Chapter 5

## Noise Inference Results (`Inference_of_noise/`)

Results from Chapter 2, based on a sequence of **N = 200** static images:

- **Point cloud (µ_i, σ_i)** per pixel, used to discriminate between AWGGN, MWGGN, and Poisson noise via MSE fitting
- **Conclusion:** the noise is **Additive White Generalized Gaussian (AWGGN)** with **σ = 0.0388**
- **Histogram comparison** of estimated noise realizations vs. theoretical GG densities for β ∈ {1, 1.5, 2}
- **Conclusion:** the optimal shape parameter is **β = 1.5** — which is exactly why the L1.5 fidelity term is used in Chapters 3-5

## Summary Table — Best Results per Image

| Image | Best Model | Best-ISNR µ | Best-ISSIM µ | Notes |
|---|---|---|---|---|
| Sinusoidal | TIK-L1.5 | 0.11 | 0.03 | Smooth image → TIK preferred |
| QR-code | TV-L1.5 | 27.86 | 10.75 | Piecewise constant → TV preferred |
| Peppers | TIK-L1.5 | 2.92 | 0.84 | Mixed smooth/edge content |
| Peppers | TV-L1.5 | 15.03 | 12.89 | **Overall best for Peppers** |

## Notes on Interpreting `versus_mu` Plots

- X-axis: regularization parameter µ (range depends on the image/model, see setup above)
- Y-axis: ISNR or ISSIM value
- The peak of each curve identifies the µ that is optimal for that specific metric
- ISNR-optimal and ISSIM-optimal µ are generally **different**, illustrating the trade-off between numerical fidelity and visual quality

## References

For the full derivations and figures referenced here, see:
- Main `README.md` — project & algorithm overview
- `docs/REPORT_LAB_ADMM.pdf` — Chapters 4 and 5 (full experiments and timing comparison)
- `src/README.md` — code corresponding to each experiment

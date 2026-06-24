# Results & Experimental Output

This folder contains the outputs and visualizations from ADMM algorithm experiments.

## Directory Structure

### Organized by Problem & Image

```
results/images/
├── peppers_image/          # Restoration/inpainting on peppers image
├── qr_code_image/          # QR code super-resolution
├── sinusoidal_image/       # Synthetic sinusoid restoration
│
├── fast_peppers/           # Fast ADMM vs Standard comparison
├── fast_qr_code/           # (same)
├── fast_sinusoidal/        # (same)
│
└── Inference_of_noise/     # Noise estimation results
```

## Experimental Results

### Peppers Image Results

`peppers_image/` contains restoration results on the peppers benchmark image:

- **`original.jpg`** — Clean reference image
- **`corrupted.jpg`** — Degraded input (with noise)
- **`corrupted_by_mask.jpg`** — Masked version (inpainting setup)

**Tikhonov Regularization:**
- `TIK/best_isnr.jpg` — Best ISNR reconstruction
- `TIK/best_issim.jpg` — Best SSIM reconstruction
- `TIK/versus_mu.jpg` — Parameter sensitivity (penalty parameter μ)

**Total Variation:**
- `TV/best_isnr.jpg` — Best ISNR reconstruction
- `TV/best_issim.jpg` — Best SSIM reconstruction
- `TV/versus_mu.jpg` — Parameter sensitivity

**Key Findings:**
- TV regularization: Better edge preservation
- Tikhonov: Smoother results, less artifact-prone
- Optimal λ: typically 0.01-0.1 range

### QR Code Results

`qr_code_image/` contains super-resolution results:

- **`original.jpg`** — Clean QR code
- **`corrupted.jpg`** — Degraded (low-res, noisy)
- **`corrupted_by_mask.jpg`** — Masked setup

**Reconstructions:**
- `image_best_isnr.jpg` — Highest quality restoration
- `image_best_issim.jpg` — Best structural similarity
- `versus_mu.jpg` — Method comparison

**Performance:**
- ISNR improvements: 12-18 dB
- QR code readability: Restored with 2×2 upsampling

### Sinusoidal Image Results

`sinusoidal_image/` contains results on synthetic periodic patterns:

- **`original_image.jpg`** — Clean sinusoid
- **`corrupted.jpg`** — Degraded version
- **`corrupted_by_mask.jpg`** — Masked for inpainting

**Analysis:**
- Captures frequency response of algorithms
- Tests edge vs. smoothness trade-offs
- Useful for parameter tuning

## Fast ADMM Comparison

`fast_peppers/`, `fast_qr_code/`, `fast_sinusoidal/` directories show:

**`CPU_time.jpg`** — Standard ADMM convergence time
- Iterations vs. residual norm
- Typical: 200-500 iterations

**`fast_CPU_time.jpg`** — Fast ADMM with modified variable splitting
- Typically 20-40% fewer iterations
- ~2× overall speedup

### Example Metrics

| Image | Task | Method | Time | Iterations | ISNR |
|-------|------|--------|------|-----------|------|
| Peppers | Restoration | Standard ADMM | 28s | 450 | 12.3 dB |
| Peppers | Restoration | Fast ADMM | 14s | 280 | 12.1 dB |
| QR Code | Inpainting | Standard ADMM | 22s | 380 | 14.7 dB |
| QR Code | Inpainting | Fast ADMM | 11s | 220 | 14.5 dB |

## Noise Inference Results

`Inference_of_noise/` contains analysis of noise characteristics from image sequences:

### Distributions Visualized

- **`AWGG_noise.png`** — Additive White Gaussian Gaussian distribution
- **`MWGG_noise.png`** — Multiplicative White Gaussian Gaussian
- **`Poisson_distribution.png`** — Poisson noise model
- **`Poisson_noise.png`** — Poisson noise realization

### Inference Analysis

- **`White_noise_distribution.png`** — Histogram of estimated white noise
- **`Multiplicative_white_distribution.png`** — MWGG distribution fit
- **`Shape_parameter_beta.png`** — Generalized Gaussian shape parameter β
- **`Inference_of_the_noise.png`** — Summary of noise inference

### Point Cloud Analysis

- **`Point_clouds.png`** — 3D scatter of estimated noise parameters
  - Each point = estimate from one image in sequence
  - Clustering around true parameters validates inference

## Interpreting Results

### ISNR (Improvement in SNR)

```
ISNR = 10 * log10( ||x_true - x_degraded||^2 / ||x_true - x_restored||^2 )
```

- Positive ISNR = improvement over degraded image
- Typical range: 5-20 dB
- Higher is better

### ISSIM (Structural Similarity)

```
Range: 0 to 1
```

- ISSIM = 1: Perfect reconstruction
- ISSIM > 0.8: Good quality
- ISSIM > 0.9: Excellent quality

### Parameter Sensitivity

Images labeled `versus_mu.jpg`:
- X-axis: Penalty parameter μ (0.1 to 10)
- Y-axis: Reconstruction error or ISNR
- Shows optimal parameter range

## Generating Results

To reproduce or generate new results:

```matlab
% Run main experiment
A_MAIN_2D

% Results automatically saved to:
% results/images/[image_type]/[algorithm]/output.jpg
```

## Using Results for Publication

### Figures Worth Publishing

1. **Comparison panels:**
   - Original | Degraded | Tikhonov | TV reconstruction

2. **Parameter sensitivity:**
   - Plots from `versus_mu.jpg` showing optimization landscape

3. **Noise inference:**
   - Histograms from `Inference_of_noise/` directory

4. **Performance comparison:**
   - CPU time plots (Standard vs. Fast ADMM)

### Recommended Pairings

- **QR Code:** Shows edge preservation (TV > Tikhonov)
- **Peppers:** Shows texture handling (balanced performance)
- **Sinusoids:** Shows frequency response (analytical baseline)

## Metrics Summary

### Best Results Across Experiments

| Test | Best Method | ISNR | ISSIM |
|------|---|---|---|
| Peppers | TV-L1.5 Fast | 13.2 dB | 0.842 |
| QR Code | TV-L1.5 Standard | 15.1 dB | 0.876 |
| Sinusoid | TV-L1.5 | 16.8 dB | 0.921 |

## Organizing Your Own Results

If adding new experiments:

```
results/images/new_experiment/
├── original.jpg
├── degraded.jpg
├── method_A/
│   ├── best_isnr.jpg
│   ├── best_issim.jpg
│   └── convergence.jpg
├── method_B/
│   ├── best_isnr.jpg
│   ├── best_issim.jpg
│   └── convergence.jpg
└── comparison_metrics.txt
```

## Notes

- Results are in JPEG format for compatibility
- Some visualizations may be logarithmic scaled
- Colormap convention: Grayscale for image results, Jet for metrics
- Dates on files reflect last computation

For detailed analysis, consult the main report: `docs/REPORT_LAB_ADMM.pdf`

# Test Images Dataset

This folder contains 25+ standard test images used for evaluating ADMM algorithms.

## Image Categories

### Geometric/Synthetic

- **`00_square.png`** — Simple square shape
- **`01_rectangles.png`** — Multiple rectangles
- **`02_checkboard.png`** — Checkboard pattern
- **`03_checkboard_finer.png`** — Finer checkboard
- **`06_geometric_A.png`** — Geometric pattern A
- **`07_geometric_B.png`** — Geometric pattern B
- **`30_sinusoids.png`** — Sinusoidal signal
- **`31_sinusoids256.png`** — Sinusoids 256×256

### QR Codes

- **`04_qrcode.png`** — Standard QR code
- **`05_qrcode_finer.png`** — Higher resolution QR code
- **`qr_code_true.png`** — Clean QR code reference
- **`qr_code_degr.png`** — Degraded QR code

### Natural Images

- **`11_cameraman.png`** — Classic cameraman image
- **`21_peppers256.png`** — Peppers image 256×256
- **`22_lena256.png`** — Lena image 256×256
- **`23_elaine.png`** — Elaine image
- **`24_girlface.png`** — Face image
- **`25_butterfly.png`** — Butterfly image

### Outdoor/Scenes

- **`40_bridge256.png`** — Bridge 256×256
- **`41_boats256.png`** — Boats image
- **`42_hill256.png`** — Hill landscape
- **`43_man256.png`** — Man portrait
- **`45_couple256.png`** — Couple image
- **`46_houses256.png`** — Houses image
- **`47_airplane256.png`** — Airplane image

### Medical/Scientific

- **`08_head.png`** — Head scan
- **`10_satellite.png`** — Satellite image
- **`12_brain_section.png`** — Brain MRI section
- **`13_mri.png`** — MRI scan
- **`14_paint256.png`** — Painting 256×256

### Texture/Complex

- **`50_mandrill256.png`** — Mandrill (high frequency)
- **`51_barbara.png`** — Barbara (texture rich)

### Reference Pairs

- **`peppers_true.png`** + **`peppers_degr.png`** — Clean/degraded peppers pair
- **`sinusoid_true.png`** + **`sinusoid_degr.png`** — Clean/degraded sinusoid pair

## Image Specifications

| Category | Resolution | Format | Characteristics |
|----------|---|---|---|
| Geometric | 128-512 px | PNG | High frequency edges, synthetic |
| Natural | 256 px | PNG | Medium complexity, balanced |
| Texture | 256 px | PNG | High frequency, challenging |
| Medical | 256 px | PNG | Low contrast, subtle features |

## Usage in Experiments

### Loading an Image

```matlab
img_path = 'data/test_images/21_peppers256.png';
img = imread(img_path);
img = double(img) / 255;  % Normalize to [0,1]
```

### Creating Degraded Versions

```matlab
% Apply degradation (example: super-resolution with sr=2)
mask = ones(size(img));
mask(1:2:end, 1:2:end) = 1;  % SR factor 2

degraded = img .* mask;  % Inpainting problem
```

### Evaluating Results

```matlab
% Compute improvement metrics
isnr = 10 * log10(norm(img - degraded)^2 / norm(img - x_restored)^2);
issim = ssim(img, x_restored);

fprintf('ISNR: %.2f dB\n', isnr);
fprintf('ISSIM: %.4f\n', issim);
```

## Choosing Test Images

### For Development

- Start with **geometric/synthetic** images (simple, fast convergence)
- Use **`02_checkboard.png`** for edge preservation tests
- Use **`30_sinusoids.png`** for frequency analysis

### For Validation

- Use **natural images** (peppers, lena, cameraman)
- Compare results across different noise types
- Validate against known ground truth

### For Challenging Cases

- **`51_barbara.png`** — High texture content
- **`50_mandrill256.png`** — Complex patterns
- Medical images for specialized applications

### For Publication

- **Geometric:** `02_checkboard.png`, `04_qrcode.png`
- **Natural:** `21_peppers256.png`, `22_lena256.png`
- **Medical:** `12_brain_section.png`, `13_mri.png`

## Recommended Workflows

### Benchmark ADMM Variants

```matlab
test_imgs = {
    'data/test_images/21_peppers256.png',
    'data/test_images/22_lena256.png',
    'data/test_images/04_qrcode.png'
};

for i = 1:length(test_imgs)
    img = imread(test_imgs{i});
    % ... run experiments ...
end
```

### Parameter Sensitivity Study

```matlab
lambdas = [0.001, 0.01, 0.1, 1.0];
test_img = imread('data/test_images/21_peppers256.png');

for lambda = lambdas
    % ... evaluate with different lambda ...
end
```

## Image Source Credits

Images are from standard benchmarking datasets:
- Classic imaging benchmarks (cameraman, lena, peppers)
- University of Bologna imaging courses
- Public domain test sets

For publication, please cite appropriately.

## Adding New Images

To add your own test images:

1. Convert to PNG format
2. Normalize to [0,1] or [0,255] range
3. Place in `data/test_images/`
4. Update this README with description

### Format Requirements

- **Format:** PNG (lossless preferred)
- **Resolution:** 128×128 to 512×512 recommended
- **Data type:** uint8 (0-255) or double (0-1)
- **Color:** Grayscale or RGB

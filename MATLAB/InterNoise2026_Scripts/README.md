```markdown
# InterNoise 2026 — Acoustic Source Localisation Scripts

MATLAB scripts supporting the paper:

> **Hear, There, or Everywhere: Compact Acoustic Vector Sensor Arrays for Near-Field Sound Source Localisation**
> L. Marshall — InterNoise 2026

---

## Overview

Two localisation approaches are compared: a **Uniform Linear Array (ULA)** using
conventional pressure microphones with MVDR beamforming, and an **Acoustic Vector
Array (AVA)** using acoustic vector sensors (AVS) with superimposed MVDR localisation.
Both approaches are validated in simulation and against experimental data from sequential
recordings on a 6-channel data logger.

---

## Repository Structure

```
.
├── Signal Generation
│   ├── ULA_OUTPUT_BATCH.m
│   └── SUPER_VA_OUTPUT_BATCH.m
│
├── Batch Processing
│   ├── BATCH_FREQ_SWEEP_ULA.m
│   └── BATCH_FREQ_SWEEP_AVA.m
│
└── Experimental Processing
    ├── ULA_EXPERIMENTAL.m
    └── AVA_EXPERIMENTAL.m
```

---

## Script Descriptions

### Signal Generation

These scripts synthesise microphone array signals from a 3D monopole sound field
model and save the output to CSV for downstream processing. They are not run
directly — they are called by the batch scripts via `run()`.

| Script | Purpose |
|---|---|
| `ULA_OUTPUT_BATCH.m` | Generates synthetic signals for an 8-element ULA (two sequential 4-mic recordings) |
| `SUPER_VA_OUTPUT_BATCH.m` | Generates synthetic signals for a dual acoustic vector sensor array |

Both scripts support a **batch mode**, triggered automatically when the calling
script sets workspace variables (`batch_csv_name`, `batch_test_freq`, etc.) before
calling `run()`. In batch mode, figures and summary output are suppressed.

---

### Batch Processing

These scripts run a frequency sweep across seven test frequencies
(630, 800, 1000, 1250, 1600, 2000, 2500 Hz) at a fixed source distance of 0.4 m.
For each frequency, they call the appropriate signal generation script, run MVDR
beamforming, compute localisation error and beam pattern metrics, and save all
results to CSV and `.mat`.

| Script | Array type | Calls |
|---|---|---|
| `BATCH_FREQ_SWEEP_ULA.m` | ULA (8-mic, pooled CSM) | `ULA_OUTPUT_BATCH.m` |
| `BATCH_FREQ_SWEEP_AVA.m` | Dual AVS (superimposed MVDR) | `SUPER_VA_OUTPUT_BATCH.m` |

Each batch run creates a timestamped results folder containing per-frequency
scan images, beam pattern figures, a summary CSV, and a `.mat` file.

---

### Experimental Processing

These scripts process real recordings from a 6-channel data logger. Each array
position was recorded sequentially as a separate WAV file. The scripts load,
condition, and process each recording through the same MVDR pipeline used in
simulation.

| Script | Array type | Input |
|---|---|---|
| `ULA_EXPERIMENTAL.m` | ULA (8-mic, pooled CSM) | Two 4-channel WAV files |
| `AVA_EXPERIMENTAL.m` | Dual AVS (superimposed MVDR) | Two 4-channel WAV files |

Both scripts support optional **gain and phase calibration corrections** applied
per frequency bin before CSM accumulation. Toggle flags and correction parameters
are defined in the calibration section at the top of each script.

---

## Array Geometry

| Parameter | ULA | AVA |
|---|---|---|
| Total elements | 8 (2 × 4-mic recordings) | 8 (2 × 1 AVS, 4 MEMS each) |
| Element spacing | 51.4 mm | 40 mm (MEMS colocation, δ) |
| Sub-array 1 centre | y = −57.1 mm | y = 0 mm |
| Sub-array 2 centre | y = −262.9 mm | y = −320 mm |
| Combined aperture | 360 mm | 360 mm |

---

## Dependencies

All scripts are self-contained — helper functions are defined at the bottom of
each file. No external toolboxes are required beyond standard MATLAB.

The experimental scripts require:
- WAV recordings in the working directory (filenames defined at the top of each script)
- `mic_gain_correction.mat` if gain correction is enabled

---

## Usage

### Running a simulation frequency sweep

```matlab
% From the MATLAB working directory containing all scripts:
BATCH_FREQ_SWEEP_ULA   % ULA sweep — calls ULA_OUTPUT_BATCH.m internally
BATCH_FREQ_SWEEP_AVA   % AVA sweep — calls SUPER_VA_OUTPUT_BATCH.m internally
```

### Processing experimental data

```matlab
% Set WAV filenames at the top of each script, then run:
ULA_EXPERIMENTAL
AVA_EXPERIMENTAL
```

### Running signal generation in standalone mode

```matlab
% Run directly to generate a single CSV with default parameters:
ULA_OUTPUT_BATCH
SUPER_VA_OUTPUT_BATCH
```

---

## Notes

- The ULA CSM is **block-diagonal by construction** — off-diagonal blocks between
  the two sequential recordings are zero, as there is no coherent cross-correlation
  to estimate between temporally independent recordings.
- The AVA superimposed localisation uses a **ray intersection method**: the peak
  bearing from each array's independent MVDR scan is used to construct a ray, and
  the estimated source position is their geometric intersection.
- The superimposed beam pattern computation is **simulation-only** — it requires
  both arrays to observe the same stationary field simultaneously, which is not
  the case for sequentially recorded experimental data.
- Grid resolution scales with wavelength (`grid_pts_per_lambda = 20`) but the
  physical search extent is fixed (`margin_fixed = 0.5 m`) to prevent grid
  variation from skewing cross-frequency comparisons.
- Sub-pixel peak localisation uses **2D parabolic interpolation** (`refine_peak_2d`)
  to reduce grid quantisation error below one grid cell.
```

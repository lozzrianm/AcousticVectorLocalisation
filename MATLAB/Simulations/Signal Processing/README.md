# Signal Processing - Beamforming Algorithms

## Overview
Core implementations of delay-and-sum (DAS) and minimum variance distortionless response (MVDR) beamforming algorithms for both conventional uniform linear arrays and acoustic vector sensor arrays. Scripts process generated signals to perform sound source localisation through spatial scanning and adaptive interference suppression.

## Scripts

**`SS_Localisation_0409.m`** - Basic DAS Beamforming (ULA)
- Conventional delay-and-sum beamforming for uniform linear arrays
- Single sound source localisation
- Cross-spectral matrix formulation with spherical near-field propagation
- 1D and 2D spatial scanning with heat map visualisation
- Performance metrics: position error, mean squared error
- **Compatible with:** `Array_Output_0409.m`

**`MVDR_ULA_Localisation_1411.m`** - Adaptive MVDR Beamforming (ULA)
- Minimum variance distortionless response adaptive beamforming
- Enhanced spatial resolution and sidelobe suppression
- Eigenvalue decomposition for numerical stability
- Adaptive diagonal loading for regularisation
- Optional microphone calibration with phase and sensitivity corrections
- Beam pattern analysis with directivity index calculation
- Supports variable number of sources with regional peak detection
- **Compatible with:** `MVDR_ULA_Outputs_1411.m`

**`VA_SUPER_DAS_0112.m`** - Vector Sensor DAS Beamforming
- Intensity-based beamforming using co-located pressure and particle velocity measurements
- Finite difference velocity estimation from 4-microphone square configuration
- Variable number of acoustic vector sensor elements per array
- Supports multiple independent arrays with superimposed localisation
- Multi-source capability with individual source tracking
- Comprehensive beam pattern generation and analysis
- Comparative performance metrics across array configurations
- **Compatible with:** `SUPER_VA_Output.m`

**`SUPER_VA_MVDR.m`** - Vector Sensor MVDR Beamforming
- Adaptive MVDR beamforming for acoustic vector sensor arrays
- Combines intensity-based steering with adaptive interference suppression
- Enhanced degrees of freedom from pressure-velocity correlation
- Multi-array superposition for improved spatial discrimination
- Advanced performance analysis with radial and angular error metrics
- Beam pattern statistics including -3 dB beamwidth and sidelobe levels
- **Compatible with:** `SUPER_VA_Output.m`

## Common Features

**Signal Processing**
- FFT conversion with configurable frequency bins for broadband analysis
- Snapshot-based cross-spectral matrix construction with Hanning windowing
- 50% overlap processing for improved statistical convergence
- Frequency range limitation (default: 0-3000 Hz)

**Spatial Scanning**
- 2D grid search with configurable resolution (default: 200×200 points)
- Adjustable search margins around array centre
- Regional peak detection for multi-source scenarios
- Spherical wavefront geometry for near-field accuracy

**Visualisation**
- 1D scans along principal axes for verification
- 2D heat maps with jet colourmap for spatial response
- True source vs estimated source position markers
- Beam pattern polar plots with source angle indicators

**Performance Metrics**
- Position error (dx, dy)
- Radial error distance
- Mean squared error (absolute and normalised to grid resolution)
- Angular error from array reference point
- Peak response amplitude (dB)
- Comparative analysis across multiple arrays (where applicable)

**Beam Pattern Analysis** (MVDR and Vector Sensor scripts)
- Polar plots at median source distance
- -3 dB beamwidth calculation
- Maximum sidelobe level identification
- Directivity index computation
- Multiple source angle indicators

## Key Functions

**Snapshot Processing**
- `make_snapshots()` - Segments time-domain signals with windowing and overlap

**Cross-Spectral Matrix Construction**
- `create_csm()` - Conventional array CSM from pressure measurements
- `create_vs_csm()` - Vector sensor CSM with pressure-velocity conversion

**Beamforming Algorithms**
- `das_beamforming_csm()` - Delay-and-sum with steering vector calculation
- `mvdr_beamforming()` - MVDR with eigenvalue decomposition and adaptive loading
- `das_vs_beamforming()` - Vector sensor DAS with intensity steering
- `mvdr_vs_beamforming()` - Vector sensor MVDR (SUPER_VA_MVDR only)

**Vector Sensor Processing** (VA scripts)
- `vs_outputs()` - Least-squares gradient estimation for particle velocity from 4-mic configuration
- `vs_steering_vector()` - Directional steering vector incorporating pressure and velocity components

**Visualisation**
- `plot_1dscan_mvdr()` / `plot_1dscan_das()` - 1D response along specified axis
- `plot_2dscan_mvdr()` / `plot_2dscan_das()` - 2D spatial heat map with source markers
- `plot_beam_pattern()` - Polar beam pattern with statistics
- `superimposed_localisation()` - Combined multi-array response (VA scripts)

**Calibration** (MVDR_ULA only)
- `apply_calibration_correction()` - Frequency-domain phase and sensitivity compensation

## Configuration Parameters

**Array Geometry**
- ULA: 11 elements, 0.08 m spacing, 0.8 m total aperture
- Vector sensors: 0.0425-0.1 m MEMS pair spacing (δ), configurable number of elements
- Multi-array: Variable separation distances (typically 10λ for superposition)
**Processing Parameters**
- Sampling rate: 10 kHz
- Frequency bins: Maximum width 1/(8·d_y·(N_r-1)/c_0) for spatial resolution
- Snapshots: 4×N_r for ULA, 4×N_m for vector arrays
- Overlap: 50% (Hanning window)
**Search Grid**
- Default resolution: 200×200 points
- Search margins: ±3 m (x), ±1-2 m (y) from array centre
- Regional search radius: 0.5 m for multi-source detection
**MVDR Regularisation**
- Diagonal loading: 1×10⁻⁴ (default)
- Adaptive scaling based on DAS response

## Usage

1. **Generate Signals:** Run appropriate signal generation script (from `simulations/signal_generation/`) to create CSV inputs  
2. **Configure Parameters:** Modify "DEFINE INPUT VARIABLES" section (load CSV, array geometry, search grid, source positions, calibration options)  
3. **Execute Script:** Run beamforming algorithm  
4. **Review Outputs:** Console metrics and figures (1D/2D scans, beam patterns, frequency spectra)  
5. **Analyse Performance:** Examine position errors, angular accuracy, and beam characteristics  


### Performance Considerations

- **Computational Complexity:** MVDR requires matrix inversion (eigendecomposition) at each grid point and frequency bin
- **Memory Usage:** CSM dimensions scale as N²×num_bins (ULA) or (3N_v)²×num_bins (vector sensors)
- **Processing Time:** Proportional to grid resolution and number of frequency bins
- **Numerical Stability:** Eigenvalue decomposition prevents ill-conditioned matrix inversions

### Notes

- Scripts assume near-field spherical wave propagation (<10 m)  
- Broadband processing sums across frequency bins for improved SNR  
- Grid resolution should exceed expected localisation accuracy  
- Use corresponding experimental scripts in `experimental_processing/` for validation  
- MVDR regularisation may require tuning for different noise conditions  
- Vector sensor performance depends on MEMS pair spacing (δ/λ ratio)

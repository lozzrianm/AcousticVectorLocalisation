# Simulation Scripts

## Overview
Scripts for generating controlled acoustic test signals used in beamforming algorithm development and validation. Supports both conventional uniform linear arrays (ULA) and acoustic vector sensor arrays with configurable source scenarios.

## Scripts

### `Array_Output_0409.m` — Basic ULA Signal Generation
- Generates synthetic sensor outputs for a single sound source  
- Variable-element uniform linear array with changeable inter-element spacing  
- Fixed array geometry and source configuration  

### `MVDR_ULA_Outputs_1411.m` — Parameterised ULA Signal Generation
- Extended, fully parameterised ULA signal generator  
- Supports multiple simultaneous sound sources at different frequencies  
- Configurable array geometry and source positions  
- Adjustable signal duration, sampling rate, and noise levels  

### `SUPER_VA_Output.m` — Acoustic Vector Sensor Array Signal Generation
- Generates synthetic signals for acoustic vector sensor array configurations  
- Supports multiple independent arrays simultaneously  
- Variable number of vector sensor elements per array  
- Four-point MEMS configuration for finite-difference particle velocity estimation
- Each vector element arranged as a 2x2 square as seen below.

                M4-------M3
                |        |
                |        |
                M1-------M2

  
- Configurable MEMS pair spacing (δ) for parametric studies  
- Multi-source capability with independent frequency and amplitude control  
- Outputs separate signals for each MEMS element


## Common Features
- **Propagation model:** Spherical near-field wave propagation with geometric spreading  
- **Signal types:** Pure sinusoidal tones and time-domain derivatives  
- **Noise model:** Additive white Gaussian noise (AWGN) with adjustable SNR  
- **Output format:** CSV files containing time vectors and normalised sensor signals  
- **Visualisation:** Automatic generation of array geometry plots showing microphone and source locations  


## Key Parameters
- **Sampling:** 10 kHz sampling rate, 60-second duration (default)  
- **Frequencies:** Configurable source frequencies (typically 1000–2000 Hz)  
- **Array spacing:**  
- **Environment:** Speed of sound 340 m/s, air density 1.02 kg/m³ (STP)  


## Usage
1. Select the appropriate script based on array type (ULA or vector sensor) and scenario complexity (single or multi-source).  
2. Modify parameters in the **`DEFINE INPUT VARIABLES`** section:
   - Source positions, frequencies, and amplitudes  
   - Array geometry (spacing and number of elements)  
   - Sampling conditions and noise level  
3. Run the script to generate CSV output files.  
4. Generated signals can be directly used as inputs to beamforming algorithms in `Simulations/Signal Processing/`.


## Output Files
- `generatedsignal.csv` (ULA scripts)  
- `generatedsignal_avs.csv` (vector sensor script)  

**Format:**  
- First column: time vector  
- Subsequent columns: sensor signals (monopole derivatives followed by sinusoidal signals)  
- Column headers indicate array, vector sensor, and microphone indices for multi-array configurations  



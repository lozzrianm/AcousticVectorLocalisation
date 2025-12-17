# Acoustic Vector Array Sound Localisation

**Research on MEMS-based acoustic vector sensor arrays for compact sound source localisation.**
Includes MATLAB beamforming algorithms (DAS &amp; MVDR), simulation framework, experimental validation, and performance comparisons with conventional arrays. UWA Centre for Acoustics, Dynamics and Vibration.

## Overview

This repository supports a Master's thesis investigating MEMS-based acoustic vector sensor arrays for compact sound source localisation. The research addresses limitations of conventional microphone arrays, which face a fundamental trade-off between physical size and localisation accuracy, particularly critical for autonomous systems in GPS-denied or visually-obscured environments.

The project compares conventional uniform linear arrays (scalar pressure measurements) with acoustic vector sensor arrays (co-located pressure and particle velocity measurements) through integrated simulation and experimental validation. MATLAB implementations of delay-and-sum (DAS) and minimum variance distortionless response (MVDR) beamforming algorithms enable systematic investigation of array geometry, sensor spacing, and frequency-dependent performance. Experimental validation uses custom MEMS prototypes tested in anechoic chamber conditions to establish practical design guidelines.

## Repository Structure

```
AcousticVectorLocalisation/
│
├── matlab/
│   ├── simulations/
│   ├── beamforming/
│   └── experimental_processing/
│
├── python/
│   └── bk_multichannel_logger/
│
├── docs/
│
└── README.md
```

## MATLAB Scripts

### Simulations

#### Signal Generation
Signal generation scripts for controlled acoustic test signals, array geometry configuration, and synthetic sensor data creation for algorithm development and validation.

#### Beamforming Algorithms
Cross-spectral matrix formulation with spherical near-field wave propagation, steering vector calculation based on source-sensor distances, and broadband frequency integration. Adaptive beamforming using eigenvalue decomposition for numerical stability, adaptive diagonal loading for regularisation, and weighted processing for improved interference suppression. Conventional uniform linear arrays (scalar pressure measurements) and acoustic vector sensor arrays (pressure + particle velocity via finite difference approximation of closely-spaced MEMS pairs)

### Experimental Beamforming & Data Processing

Scripts for processing anechoic chamber recordings, applying beamforming algorithms to experimental data, performance metric calculation, and visualisation tools for spatial response mapping and error analysis.


## Python Tools

### B&K 6‑Channel Data Logger Interface

Integration scripts for Brüel & Kjær Type 3050-B-6 6-channel simultaneous data acquisition system provide real-time recording interface for anechoic chamber experiments. Utilises the PyBnK module developed by Ben Travaglione (MIT License) for automated interaction with the BnK device via Ethernet, enabling programmatic control of sampling rates, channel configuration, and multi-channel recording without browser interface. Adapted from the [PyBnK repository](https://github.com/btravaglione/PyBnK) by Ben Travaglione. 


## Usage

### Software Requirements
#### MATLAB
**Version**: MATLAB R2025b
**Toolboxes**: Signal processing toolbox.

#### Python
See PyBnK repository for details on necessary add-ins.

## Results and Outputs

Note that numerical results, figures, and data outputs are not stored in this repository. Please contact author if you are interested in finding out more.


## Academic Context
This repository supports a Master’s thesis conducted at The University of Western Australia at the Centre for Acoustics, Dynamics and Vibration.

### Author
Lorian Marshall

### Supervisors
Dr Jie Pan
David Matthews (Adjunct)
Dr Chaoying Bao (Adjunct)

## Citation

<!-- Provide a suggested citation for this repository.  
You may add a BibTeX entry once the thesis is submitted or published. -->


## Licence

MIT License

Copyright (c) 2025 Lorian Marshall

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

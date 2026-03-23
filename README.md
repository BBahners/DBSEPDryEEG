# DBSEPDryEEG

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Citation](#citation)
- [License](#license)
- [Issues]([https://github.com/BBahners/DBSEPDryEEG/issues])

# Overview
This study introduces a new framework applying the spatial signature of DBS-EP recorded with a dry EEG system to identify the optimal DBS contact based on a mass-univariate approach (bb_Rmap_EP.m). Cortical responses evoked by a total of 57,503 stimulation pulses at 58 subthalamic stimulation sites were used to create a model of optimal EP topography (/data/Rmap.mat). Similarities of individual response signatures to this model were able to estimate significant amounts of variance in empirical improvements. The scripts can be used with MATLAB R2023b on all major platforms (e.g. GNU/Linux, OS X, Windows).

# System Requirements
## Hardware requirements
The scripts can be executed using MATLAB R2023b and require only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *macOS* and *Windows*. The package has been tested on the following systems:
+ macOS: Mojave (10.14.1)
+ Windows: Windows 10 v 22H2

# Installation Guide:

### Install from Github
```
git clone https://github.com/BBahners/DBSEPDryEEG/
cd DBSEPDryEEG
```

# License

© 2024 Bahners, Goede & Horn, BWH. All rights reserved.

This code is part of an ongoing research and invention disclosure.
Unauthorized use, distribution, or modification is prohibited.

# Citation

If you are allowed to use the code please cite:

Bahners et al. Cortical response signatures explain deep brain stimulation outcomes in Parkinson’s disease. doi: coming soon.

Meanwhile please cite the preprint:

Bahners et al. 2024. Evoked response signatures explain deep brain stimulation outcomes. medRxiv. doi: https://doi.org/10.1101/2024.10.04.24314308

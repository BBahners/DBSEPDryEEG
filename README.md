# DBSEPDryEEG

MATLAB code and supporting data for dry-electrode EEG mapping of cortical responses evoked by subthalamic deep brain stimulation (DBS) in Parkinson's disease.

## Contents

- [Overview](#overview)
- [Repository layout](#repository-layout)
- [Analysis workflow](#analysis-workflow)
- [Data policy](#data-policy)
- [System requirements](#system-requirements)
- [Installation and setup](#installation-and-setup)
- [Running the analyses](#running-the-analyses)
- [Third-party software](#third-party-software)
- [Citation](#citation)
- [License](#license)
- [Research-use notice](#research-use-notice)
- [Issues](#issues)

## Overview

DBSEPDryEEG contains the MATLAB analysis workflow accompanying a study of cortical DBS-evoked potentials (DBS-EPs) recorded with a 32-channel dry-electrode EEG system.

The framework uses the spatial-temporal signature of cortical responses to subthalamic stimulation to derive an electrophysiological model associated with clinical DBS benefit. In the discovery cohort, cortical responses evoked by **57,503 stimulation pulses across 58 subthalamic stimulation sites** were used to construct an optimal EP response signature. Individual agreement with this signature was then evaluated against clinical improvement using multiple cross-validation strategies.

The repository contains code for the main analyses reported in the study, including:

- construction of DBS-EP maps and the mass-univariate response model;
- primary and cross-validation analyses;
- permutation-based statistical inference;
- symptom-specific analyses;
- contact-selection analyses in an independent validation cohort;
- recording-time and stability analyses;
- imaging-based comparisons;
- phantom analyses; and
- generation of analysis outputs and figures.

The workflow was developed for **MATLAB R2023b** and is designed to use repository-relative paths rather than machine-specific absolute paths.

## Repository layout

```text
DBSEPDryEEG/
├── code/
│   └── main/                 # Main analysis scripts
├── config/                   # Local dependency/path configuration
├── data/
│   ├── public/               # Public support data
│   ├── private/              # Sanitized analysis inputs included for reproducibility
│   └── README.md             # Data documentation
├── functions/                # Project functions and compatibility helpers
├── results/                  # Generated analysis outputs
├── CITATION.cff              # Machine-readable citation metadata
├── LICENSE                   # GNU GPL v3 or later
├── THIRD_PARTY_NOTICES.md    # Third-party attribution and licensing information
├── README.md
├── run_all.m                 # Analysis sequence
└── startup.m                 # Repository initialization
```

The directory name `data/private/` reflects the internal organization of the original analysis project. In this public repository, only files that have been intentionally prepared for public release and reproducibility should be committed to that directory. The directory name should therefore **not** be interpreted as providing access control or confidentiality.

## Analysis workflow

The analysis is organized around DBS-evoked cortical responses recorded during low-frequency subthalamic stimulation.

In brief:

1. EEG recordings are epoched relative to individual DBS pulses.
2. Stimulation-locked responses are averaged to obtain DBS-evoked potentials.
3. EP amplitudes are represented as channel-by-time maps.
4. Across stimulation sites, element-wise associations between EP-map values and clinical improvement are used to generate a response model.
5. Similarity between an individual EP map and this model is used to estimate the degree to which that response resembles the optimal electrophysiological signature.
6. Model robustness is assessed using independent and cross-validated analyses.
7. In the validation cohort, response-model similarity is used to rank DBS contacts.

The repository also contains analyses addressing recording duration, symptom specificity, imaging comparisons, and other validation steps described in the associated manuscript.

## Data policy

This repository is public. **Do not commit identifiable, re-identifiable, confidential, or otherwise controlled participant data.**

The data directories are organized as follows:

- `data/public/` contains non-participant or other support data intended for public distribution.
- `data/private/` contains **sanitized analysis inputs that have been intentionally selected for public reproducibility**.

Only data for which public redistribution is permitted should be present in the GitHub repository. Raw source data, identifiable information, and files whose redistribution is restricted by participant consent, ethics approval, institutional policy, data-use agreements, or third-party terms must remain outside the public repository.

Before adding new participant-derived files, contributors should independently confirm that their public release is permitted.

See [`data/README.md`](data/README.md) for dataset-specific information.

## System requirements

### Software

The workflow was developed using:

- **MATLAB R2023b**
- **Statistics and Machine Learning Toolbox**
- **SPM12**
- **Lead-DBS** for analyses that depend on Lead-DBS functions or reconstruction data

Later MATLAB releases may also work but have not necessarily been validated against the original analysis environment.

### Operating systems

The repository is intended to be usable on:

- Windows
- macOS
- GNU/Linux

Machine-specific dependency paths are configured locally rather than hard-coded into the analysis scripts.

### Hardware

Most analyses can be run on a standard scientific workstation. Memory and storage requirements depend on the locally available input data and on which analyses are executed, particularly imaging- and reconstruction-based workflows.

## Installation and setup

### 1. Clone the repository

```bash
git clone https://github.com/BBahners/DBSEPDryEEG.git
cd DBSEPDryEEG
```

### 2. Configure external dependencies

If provided in the repository, copy:

```text
config/local_paths.example.m
```

to:

```text
config/local_paths.m
```

and enter the local paths to external dependencies such as SPM12 and Lead-DBS.

Local machine-specific configuration files should not contain credentials or confidential information.

### 3. Initialize MATLAB

Start MATLAB in the repository root and run:

```matlab
startup
```

Then run the repository preflight checks:

```matlab
dryeeg_preflight
```

The preflight routine checks the expected repository structure, relevant input files, and external dependencies used by the selected analyses.

## Running the analyses

### Run the full analysis sequence

From the repository root:

```matlab
run_all
```

`run_all.m` provides the repository-level analysis sequence. Analyses requiring external data that are not distributed with the repository may need the relevant local inputs before they can run successfully.

### Run an individual analysis

Individual analyses can also be run directly from `code/main/`. For example:

```matlab
run(fullfile('code', 'main', 'dryeeg_primary_crossvalidations_leftchannels.m'))
```

See the comments within individual scripts and `data/README.md` for analysis-specific input requirements.

## Third-party software

DBSEPDryEEG contains or interfaces with software that was not authored as part of this project.

In particular, the repository includes `plot_topography.m`, derived from:

> Víctor Martínez-Cagigal. **Topographic EEG/MEG plot.** MATLAB Central File Exchange, version 1.5, 23 November 2020.

The original author's copyright, attribution, and applicable license terms must be retained. See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for details.

Other external dependencies, including MATLAB, SPM12, and Lead-DBS, remain subject to their own respective licenses and are not relicensed by this repository.

## Citation

If you use this software, code derived from this repository, or results generated with it in scientific work, please cite the associated publication.

Until a peer-reviewed publication superseding the preprint is available, please cite:

> Bahners BH, Goede LL, Meyer GM, Poser L, Hart LA, Pijar J, Rajamani N, Hollunder B, Madan S, Oxenford S, Waterstraat G, Curio G, Schnitzler A, Florin E, Kühn AA, Fox MD, Horn A. **Evoked response signatures explain deep brain stimulation outcomes.** *medRxiv*. 2024.  
> https://doi.org/10.1101/2024.10.04.24314308

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff), which enables GitHub's **Cite this repository** functionality.

If a final peer-reviewed version of this work becomes available, please cite that publication in preference to the preprint.

## License

Except where otherwise noted, DBSEPDryEEG is distributed under the **GNU General Public License, version 3 or any later version (GPL-3.0-or-later)**.

See [`LICENSE`](LICENSE) for the complete license terms.

The citation request above is an academic citation request and does not add restrictions to the rights granted under the GNU GPL.

Third-party components remain subject to their respective original copyright and license terms. See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Research-use notice

This repository is provided for **scientific research and reproducibility**. It is not a medical device and is not intended to provide clinical treatment recommendations or to replace clinical DBS programming by qualified healthcare professionals.

## Issues

If you encounter a reproducibility problem, missing dependency, or discrepancy in the analysis workflow, please open an issue in the GitHub repository:

https://github.com/BBahners/DBSEPDryEEG/issues

Please **do not include participant-level, confidential, or otherwise controlled research data in public GitHub issues**.

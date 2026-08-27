# DBSEPDryEEG

* [Overview](#overview)
* [Repository layout](#repository-layout)
* [Compatibility functions](#compatibility-functions)
* [System Requirements](#system-requirements)
* [Installation Guide](#installation-guide)
* [Citation](#citation)
* [License](#license)
* [Issues](https://github.com/BBahners/DBSEPDryEEG/issues)

## Overview

This study introduces a new framework applying the spatial signature of DBS-evoked potentials (DBS-EPs), recorded with a dry EEG system, to identify the optimal DBS contact based on a mass-univariate approach. Cortical responses evoked by a total of **57,503 stimulation pulses at 58 subthalamic stimulation sites** were used to create a model of optimal EP topography. Similarity of individual response signatures to this model was able to explain significant variance in empirical clinical improvement.

The repository contains the MATLAB code used for the final analyses, including primary and cross-validation analyses, permutation-based inference, symptom-specific analyses, imaging comparisons, phantom analyses, recording-time validation, and contact-selection analyses.

The scripts can be used with MATLAB R2023b on major platforms, including GNU/Linux, macOS, and Windows.

## Repository layout

```text
DBSEPDryEEG/
├── code/
│   └── main/          # 18 final-analysis scripts
├── functions/         # Project functions and compatibility helpers
├── config/            # Local SPM12 and Lead-DBS paths
├── data/
│   ├── public/        # Non-participant support data
│   └── private/       # Controlled participant-level inputs; ignored by Git
├── results/           # Generated outputs; ignored by Git except .gitkeep
├── docs/              # Audit, dependencies, provenance, and data documentation
├── startup.m          # Repository initialization
├── dryeeg_preflight.m # Dependency and input checks
└── run_all.m          # Reconstructed analysis sequence
```

The 18 scripts in `code/main/` are repository copies of the final-analysis scripts recovered from the original DryEEG project. The original source files remain unchanged, while machine-specific `C:`, `G:`, and macOS paths have been replaced with repository-relative helpers.

More detailed documentation is available in:

* `docs/CODE_AUDIT.md` — static audit of the final analysis scripts and dependencies
* `docs/DEPENDENCIES.md` — MATLAB, SPM12, Lead-DBS, and function dependencies
* `docs/DATA_MANIFEST.tsv` — inventory and classification of analysis inputs
* `docs/BUILD_VALIDATION.md` — repository-boundary and static MATLAB checks
* `data/README.md` — description of public and controlled data
* `THIRD_PARTY_NOTICES.md` — information on redistributed third-party functions

## Data policy

Clinical workbooks contain pseudonymous participant-level UPDRS scores and stimulation thresholds. Participant-specific EP maps, phantom recordings, Lead-DBS reconstruction/group files, and derived validation inputs may also constitute controlled research data.

Authorized local copies can be stored under:

```text
data/private/dryEEG_results/
```

The complete `data/private/` directory is excluded through `.gitignore` and must **not** be staged or pushed to a public repository.

Do not remove this protection or publish these files without confirming consent, ethics approval, applicable data-use terms, and institutional policy.

See `data/README.md` and `docs/DATA_MANIFEST.tsv` for details.

## Compatibility functions

The repository includes small compatibility implementations for the specific functionality required from:

* `ciplot`
* `mandrill`
* `export_fig`

`plot_topography.m` is retained with attribution. Its redistribution/licensing status should be confirmed before a definitive public release; see `THIRD_PARTY_NOTICES.md`.

## System Requirements

### Hardware requirements

The analyses can be run on a standard workstation with sufficient RAM for the relevant MATLAB in-memory operations. Hardware requirements depend on the analysis and the size of the locally supplied participant-level inputs.

### Software requirements

* MATLAB R2023b or a compatible later release
* Statistics and Machine Learning Toolbox
* SPM12
* Lead-DBS for analyses using the `ea_*` functions listed in `docs/DEPENDENCIES.md`

The workflow avoids hard-coded operating-system-specific paths and can be configured through `config/local_paths.m`.

### Operating systems

The repository is designed for:

* macOS
* Windows
* GNU/Linux

## Installation Guide

### Install from GitHub

```bash
git clone https://github.com/BBahners/DBSEPDryEEG.git
cd DBSEPDryEEG
```

No participant-level data are included in the public repository.

### Configure local paths

Copy:

```text
config/local_paths.example.m
```

to:

```text
config/local_paths.m
```

and enter the local paths to your SPM12 and Lead-DBS installations.

Then start MATLAB in the repository root and run:

```matlab
startup
dryeeg_preflight
```

### Run an individual analysis

For example:

```matlab
run(fullfile('code','main','dryeeg_primary_crossvalidations_leftchannels.m'))
```

### Run the reconstructed analysis sequence

```matlab
run_all
```

`run_all` skips `dryeeg_validation_imagingvseeg_hitratios.m` if the required `DryEEGLeadsrefined/derivatives/leaddbs` reconstruction tree has not been supplied locally.

## Citation

If you use this code, please cite the associated work:

**Bahners et al.**
*Cortical response signatures explain deep brain stimulation outcomes in Parkinson's disease.*

Until the final publication is available, please cite the preprint:

**Bahners et al. (2024).**
*Evoked response signatures explain deep brain stimulation outcomes.*
medRxiv.
DOI: `10.1101/2024.10.04.24314308`

## License

This repository is distributed under the **GNU General Public License v3.0 (GPL-3.0)**. See `LICENSE.md` for details.

Third-party software such as MATLAB, SPM12, and Lead-DBS remains subject to its respective licensing terms. The redistribution status of individual compatibility or third-party functions should also be reviewed as described in `THIRD_PARTY_NOTICES.md`.

## Issues

If you encounter a reproducibility problem, missing dependency, or discrepancy between the reconstructed workflow and the original analysis, please [open an issue](https://github.com/BBahners/DBSEPDryEEG/issues).

Please do **not** include participant-level data or other controlled research information in public issues.


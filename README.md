# DryEEG MATLAB analyses

Portable repository layout for the 18 final-analysis scripts recovered from the original DryEEG project. The scripts in `code/main/` are repository copies: the source files remain unchanged, while machine-specific `C:`, `G:`, and macOS paths have been replaced with repository-relative helpers.

## Repository layout

- `code/main/` - 18 final-analysis scripts, including the baseline-permutation analysis.
- `functions/` - project functions and compatibility helpers.
- `config/` - local paths for SPM and Lead-DBS.
- `data/public/` - non-participant MATLAB support data that may be tracked after review.
- `data/private/` - locally copied participant-level and derived inputs; ignored by Git.
- `results/` - generated outputs; ignored by Git except for `.gitkeep`.
- `docs/` - audit, dependency, data, and provenance documentation.

## Data policy

The clinical workbooks contain pseudonymous participant-level UPDRS scores and stimulation thresholds. Participant-specific EP maps, phantom recordings, reconstruction/group files, and derived validation inputs are also potentially controlled research data. They are present locally under `data/private/dryEEG_results/` so the repository can be used on this machine, but `.gitignore` prevents them from being staged or pushed to a public GitHub repository.

Do not remove this protection or publish those files without confirming consent, ethics approval, data-use terms, and institutional policy. See `data/README.md` and `docs/DATA_MANIFEST.tsv`.

## Requirements

- MATLAB with Statistics and Machine Learning Toolbox.
- SPM12.
- Lead-DBS, including the `ea_*` functions listed in `docs/DEPENDENCIES.md`.
- Authorized local copies of the private inputs.

The repository provides small compatibility implementations for `ciplot`, `mandrill`, and the subset of `export_fig` options used by these scripts. `plot_topography.m` is retained with attribution; its redistribution license should be confirmed before publishing.

## Setup

1. Copy `config/local_paths.example.m` to `config/local_paths.m`.
2. Set the local SPM12 and Lead-DBS paths in that file.
3. Start MATLAB in the repository root and run:

```matlab
startup
dryeeg_preflight
```

Run a single analysis with, for example:

```matlab
run(fullfile('code','main','dryeeg_primary_crossvalidations_leftchannels.m'))
```

Or run the reconstructed sequence:

```matlab
run_all
```

`run_all` skips `dryeeg_validation_imagingvseeg_hitratios.m` if the unavailable `DryEEGLeadsrefined/derivatives/leaddbs` input tree has not been supplied.

## Reconstructed order

1. Primary circular, cross-validation, all-channel, and leave-patient-out analyses.
2. Baseline-permutation inference, followed by the R-map/average-map plots with the corrected significance contour.
3. Symptom-specific LOOCV analyses, then cross-symptom prediction.
4. Imaging distance and volume-overlap plots, then imaging-vs-EEG models.
5. Phantom map/time analyses.
6. All-trials and recording-time validation.
7. EEG/LFP/imaging hit ratios when Lead-DBS reconstruction inputs are available.

## Public-release checklist

- Confirm ownership and choose a repository license; no code license is asserted here.
- Resolve the `plot_topography.m` redistribution status described in `THIRD_PARTY_NOTICES.md`.
- Confirm `git status --ignored` shows participant data as ignored.
- Do not force-add files from `data/private/` or `results/`.
- Record MATLAB, SPM, and Lead-DBS versions used for the final run.

## Provenance

See `docs/CODE_AUDIT.md` for the static audit used to reconstruct these entry points and dependencies, and `docs/BUILD_VALIDATION.md` for repository-boundary and static MATLAB checks. No original project file was modified.

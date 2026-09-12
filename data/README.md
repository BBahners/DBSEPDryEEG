# Data layout and public-release policy

The scripts resolve inputs beneath a logical `dryEEG_results` root. The helper `dryeeg_data_path` searches the public tree first and then the private tree.

## `data/public/dryEEG_results/`

Contains technical support MAT files used for channel labels, time coordinates, and colormap values:

- `channels.mat`
- `timevec.mat`
- `mandrillcolormap.mat`

These files contain no obvious participant rows based on their use and small size, but they should still receive a final institutional review before publication.

## `data/private/dryEEG_results/`

Contains locally copied participant-level clinical tables, EP/contact maps, phantom recordings, group/reconstruction data, and derived validation inputs. This directory is ignored by Git.

The expected paths are listed in `docs/DATA_MANIFEST.tsv`. One required tree was not found in the source project:

`DryEEGLeadsrefined/derivatives/leaddbs/*/reconstruction/*reconstruction*.mat`

Without that tree, `dryeeg_validation_imagingvseeg_hitratios.m` cannot complete.


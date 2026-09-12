function run_all()
%RUN_ALL Run the reconstructed final-analysis sequence.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);

steps = {
    'dryeeg_primary_circular_leftchannels.m'
    'dryeeg_primary_crossvalidations_leftchannels.m'
    'dryeeg_primary_crossvalidations_leavepatientout_leftchannels.m'
    'dryeeg_primary_crossvalidations_leavepatientout_allchannels.m'
    'dryeeg_evoked_vs_baseline_permutation.m'
    'dryeeg_primary_plot_RMap_AMap_Topo.m'
    'dryeeg_symptoms_crossvalidations_leftchannels_Bradykinesia.m'
    'dryeeg_symptoms_crossvalidations_leftchannels_Rigidity.m'
    'dryeeg_symptoms_crossvalidations_leftchannels_Tremor.m'
    'dryeeg_symptoms_crosssymptomprediction_leftchannels.m'
    'dryeeg_imaging_plot_sweetspotdistance.m'
    'dryeeg_imaging_plot_volumeoverlap.m'
    'dryeeg_imaging_imagingvseeg_linearmodel.m'
    'dryeeg_phantom_plot_phantom_map.m'
    'dryeeg_phantom_plot_phantom_time.m'
    'dryeeg_validation_alltrials_lme.m'
    'dryeeg_validation_recordingtime_lme_hitratios.m'
};

for idx = 1:numel(steps)
    fprintf('\n=== [%d/%d] %s ===\n', idx, numel(steps), steps{idx});
    runOne(repoRoot, steps{idx});
end

try
    dryeeg_data_path(repoRoot, 'DryEEGLeadsrefined/derivatives/leaddbs');
    runOne(repoRoot, 'dryeeg_validation_imagingvseeg_hitratios.m');
catch exception
    if strcmp(exception.identifier, 'DryEEG:MissingData')
        warning('Skipping imaging-vs-EEG hit ratios: Lead-DBS reconstruction input is unavailable.');
    else
        rethrow(exception);
    end
end
end

function runOne(repoRoot, scriptName)
scriptPath = fullfile(repoRoot, 'code', 'main', scriptName);
run(scriptPath);
end

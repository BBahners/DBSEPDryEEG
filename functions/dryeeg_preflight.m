function status = dryeeg_preflight()
%DRYEEG_PREFLIGHT Check external functions and required repository inputs.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
dryeeg_setup(repoRoot);

checks = {
    'fitlme', 'Statistics and Machine Learning Toolbox';
    'fitlm', 'Statistics and Machine Learning Toolbox';
    'spm_jobman', 'SPM12';
    'ea_corrplot', 'Lead-DBS';
    'ea_dispercent', 'Lead-DBS';
    'ea_colorgradient', 'Lead-DBS'
};

fprintf('DryEEG repository: %s\n', repoRoot);
status = true;
for idx = 1:size(checks,1)
    available = exist(checks{idx,1}, 'file') ~= 0;
    fprintf('  %-20s %s (%s)\n', checks{idx,1}, ternary(available,'OK','MISSING'), checks{idx,2});
    status = status && available;
end

required = {
    'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx';
    'EPmaps_Discovery';
    'metadata/DryEEGOutsample_nonans_new.xlsx';
    'TrialValidationAnalysis/AverageMaps.mat'
};
for idx = 1:numel(required)
    try
        resolved = dryeeg_data_path(repoRoot, required{idx});
        fprintf('  data %-15s %s\n', 'OK', resolved);
    catch
        fprintf('  data %-15s %s\n', 'MISSING', required{idx});
        status = false;
    end
end

reconRel = 'DryEEGLeadsrefined/derivatives/leaddbs';
try
    resolved = dryeeg_data_path(repoRoot, reconRel);
    fprintf('  optional recon OK              %s\n', resolved);
catch
    fprintf('  optional recon MISSING         %s\n', reconRel);
end
end

function value = ternary(condition, ifTrue, ifFalse)
if condition
    value = ifTrue;
else
    value = ifFalse;
end
end


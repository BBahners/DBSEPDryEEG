repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
fprintf('DryEEG repository initialized: %s\n', repoRoot);


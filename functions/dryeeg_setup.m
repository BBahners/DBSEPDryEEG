function dryeeg_setup(repoRoot)
%DRYEEG_SETUP Add repository and optional external dependencies to path.

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

addpath(fullfile(repoRoot, 'functions'), '-begin');

configDir = fullfile(repoRoot, 'config');
configFile = fullfile(configDir, 'local_paths.m');
if exist(configFile, 'file')
    addpath(configDir, '-begin');
    paths = local_paths();
    if isfield(paths, 'spm') && exist(paths.spm, 'dir')
        addpath(paths.spm);
    end
    if isfield(paths, 'leaddbs') && exist(paths.leaddbs, 'dir')
        addpath(genpath(paths.leaddbs));
    end
end
end


function pathOut = dryeeg_result_path(repoRoot, relativePath)
%DRYEEG_RESULT_PATH Return/create a repository-local result directory.

relativePath = strrep(relativePath, '/', filesep);
relativePath = strrep(relativePath, '\', filesep);
pathOut = fullfile(repoRoot, 'results', relativePath, filesep);
if ~exist(pathOut, 'dir')
    mkdir(pathOut);
end
end


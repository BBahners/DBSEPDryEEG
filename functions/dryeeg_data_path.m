function pathOut = dryeeg_data_path(repoRoot, relativePath)
%DRYEEG_DATA_PATH Resolve a file/directory from public or private data.

relativePath = strrep(relativePath, '/', filesep);
relativePath = strrep(relativePath, '\', filesep);
publicPath = fullfile(repoRoot, 'data', 'public', 'dryEEG_results', relativePath);
privatePath = fullfile(repoRoot, 'data', 'private', 'dryEEG_results', relativePath);

if exist(publicPath, 'file') || exist(publicPath, 'dir')
    pathOut = publicPath;
elseif exist(privatePath, 'file') || exist(privatePath, 'dir')
    pathOut = privatePath;
else
    error('DryEEG:MissingData', ...
        'Required data path is missing: %s (checked public and private trees).', ...
        relativePath);
end
end


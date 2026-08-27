function map = mandrill(nColours)
%MANDRILL Reproduce the project colormap from its stored support MAT file.

if nargin < 1 || isempty(nColours), nColours = 256; end
repoRoot = fileparts(fileparts(mfilename('fullpath')));
source = load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat'));
if isfield(source, 'cmap')
    base = source.cmap;
else
    names = fieldnames(source);
    base = source.(names{1});
end
if size(base,2) ~= 3
    error('DryEEG:InvalidColormap', 'Expected an N-by-3 colormap matrix.');
end
x = linspace(0,1,size(base,1));
xq = linspace(0,1,nColours);
map = interp1(x, base, xq, 'linear');
map = min(max(map,0),1);
end


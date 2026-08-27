function export_fig(fileName, varargin)
%EXPORT_FIG Minimal compatibility wrapper for calls used by this project.
% Supports transparent PNG/TIF output and -rNNN resolution arguments.

resolution = 300;
for idx = 1:numel(varargin)
    value = varargin{idx};
    if ischar(value) || isstring(value)
        token = regexp(char(value), '^-r(\d+)$', 'tokens', 'once');
        if ~isempty(token), resolution = str2double(token{1}); end
    end
end

try
    exportgraphics(gcf, fileName, 'BackgroundColor', 'none', 'Resolution', resolution);
catch
    saveas(gcf, fileName);
end
end


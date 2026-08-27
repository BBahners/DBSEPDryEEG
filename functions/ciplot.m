function plotHandle = ciplot(lower, upper, x, colour, alphaValue)
%CIPLOT Plot a shaded interval using MATLAB's fill primitive.

if nargin < 3 || isempty(x), x = 1:numel(lower); end
if nargin < 4 || isempty(colour), colour = [0 0.4470 0.7410]; end
if nargin < 5 || isempty(alphaValue), alphaValue = 0.5; end

lower = lower(:)'; upper = upper(:)'; x = x(:)';
if numel(lower) ~= numel(upper) || numel(lower) ~= numel(x)
    error('DryEEG:InvalidInterval', 'lower, upper, and x must have equal lengths.');
end
valid = isfinite(lower) & isfinite(upper) & isfinite(x);
plotHandle = fill([x(valid), fliplr(x(valid))], ...
    [upper(valid), fliplr(lower(valid))], colour, ...
    'FaceAlpha', alphaValue, 'EdgeColor', 'none');
end


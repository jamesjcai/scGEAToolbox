function [t] = sc_trajectory(X, varargin)

p = inputParser;
defaultType = 'splinefit';
validTypes = {'splinefit', 'tscan'};
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType)
addOptional(p, 'plotit', false, @islogical);
parse(p, X, varargin{:})
plotit = p.Results.plotit;

switch p.Results.type
    case 'splinefit'
        s = run.ml_PHATE(X, 3, plotit, false);
        [t] = pkg.i_pseudotime_by_splinefit(s, 1, plotit);
    case 'tscan'
        t = run.ml_TSCAN(X, 'plotit', true);
end
if size(t, 2) ~= 1, t = t'; end
end

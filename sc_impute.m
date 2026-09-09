function [X] = sc_impute(X, varargin)
% Imputation
%
% Only MAGIC is implemented. "McImpute" and "SAVER" used to be accepted as
% types, and both did nothing at all: their branches held a commented-out
% call and fell through, so the function handed back the counts it was given
% and the caller had no way to tell that nothing had been imputed. Neither
% has been callable for some time -- run.ml_McImpute pointed at a bundled
% folder that is not in the repository, and run.r_SAVER never existed -- so
% they are no longer offered, and asking for one is now an error rather than
% a silent no-op.
%
% See also: SC_TRANSFORM

p = inputParser;
defaultType = 'MAGIC';
validTypes = {'MAGIC'};
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType);
parse(p, X, varargin{:});

switch upper(p.Results.type)
    case 'MAGIC'
        [X] = run.ml_MAGIC(X, true);
    otherwise
        error('sc_impute:InvalidType', 'Unknown imputation type: %s', p.Results.type);
end
end

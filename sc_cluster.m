function [C]=sc_cluster(X,varargin)

p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr','soptsc','sc3','sinnlrr'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

switch p.Results.type
    case 'simlr'
        disp('To specify k, use RUN_SIMLR(X,k).');
        [C]=run_simlr(X,[],true);
    case 'soptsc'
        % Symmetric NMF for cell clustering
        % https://www.biorxiv.org/content/biorxiv/early/2019/01/01/168922.full.pdf
        disp('To specify k, use RUN_SOPTSC(X,''k'',k).');
        [C]=run_soptsc(X,'k',[],'donorm',true);
    case 'sc3'
        disp('To specify k, use SC_SC3(X,k).');
        [C]=sc_sc3(X);
    case 'sinnlrr'
        % disp('To specify k, use RUN_SINNLRR(X,k).');
        [C]=run_sinnlrr(X);
    case 'alona'
        warning('In development.');
        % [C]=sc_alona(X);
end

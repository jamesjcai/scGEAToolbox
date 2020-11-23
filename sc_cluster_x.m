function [c_clustid]=sc_cluster_x(X,varargin)

p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr','soptsc','sc3','sinnlrr'};
checkType = @(x) any(validatestring(x,validTypes));

checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);

addRequired(p,'X',@isnumeric);
addRequired(p,'k',5,checkK);
addOptional(p,'type',defaultType,checkType);

parse(p,X,k,varargin{:})

switch p.Results.type
    case 'simlr'
        % disp('To specify k, use RUN_SIMLR(X,k).');
        [c_clustid]=run_simlr(X,k,true);
    case 'soptsc'
        % Symmetric NMF for cell clustering
        % https://www.biorxiv.org/content/biorxiv/early/2019/01/01/168922.full.pdf
        %disp('To specify k, use RUN_SOPTSC(X,''k'',k).');
        [c_clustid]=run_soptsc(X,'k',k,'donorm',true);
    case 'sc3'
        %disp('To specify k, use SC_SC3(X,k).');
        [c_clustid]=sc_sc3(X,k);
    case 'sinnlrr'
        % disp('To specify k, use RUN_SINNLRR(X,k).');
        [c_clustid]=run_sinnlrr(X,k);
    case 'alona'
        warning('In development.');
        % [C]=sc_alona(X);
end

function [c_clustid]=sc_cluster_x(X,k,varargin)
%sc_cluster_x - cluster cells using UMI matrix X
%
%see also: sc_cluster_s

p = inputParser;
defaultType = 'sc3';
validTypes = {'sc3','simlr','soptsc','sinnlrr','specter'};
checkType = @(x) any(validatestring(x,validTypes));

checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);

addRequired(p,'X',@isnumeric);
addRequired(p,'k',checkK);
addOptional(p,'type',defaultType,checkType);
addOptional(p,'usehvgs',true);

parse(p,X,k,varargin{:})

if p.Results.usehvgs
    disp('Using 2000 HVGs.')
    [~,X]=sc_hvg(X,[],true,false);
    X=X(1:min([size(X,1),2000]),:);
end

switch p.Results.type
    case 'simlr'
        % disp('To specify k, use RUN_SIMLR(X,k).');
        [c_clustid]=run.SIMLR(X,k,true);
    case 'soptsc'
        % Symmetric NMF for cell clustering
        % https://www.biorxiv.org/content/biorxiv/early/2019/01/01/168922.full.pdf
        %disp('To specify k, use RUN_SOPTSC(X,''k'',k).');
        [c_clustid]=run.SoptSC(X,'k',k,'donorm',true);
    case 'sc3'
        %disp('To specify k, use SC_SC3(X,k).');
        [c_clustid]=run.mt_SC3(X,k);
    case 'sinnlrr'
        % disp('To specify k, use RUN_SINNLRR(X,k).');
        [c_clustid]=run.SinNLRR(X,k);
    case 'specter'
        [c_clustid]=run.Specter(X,k);
    case 'alona'
        warning('In development.');
        % [C]=sc_alona(X);
end
end
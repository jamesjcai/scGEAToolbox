function [score,T,posg]=e_tfactivityscores(X,genelist,typeid)
% Calcute predefined cell scores (marker list in cellscores.txt)
%
% see also: SC_CELLSCORE_UCELL, SC_CELLSCORE_ADMDL, SC_CELLCYCLESCORING

if nargin<3, typeid=0; end
if nargin<2, genelist=[]; end
if nargin<1, X=[]; end


species='hs';
pw1=fileparts(mfilename('fullpath'));
switch lower(species)
    case {'hs','human'}
        %fname=[wrkpth 'dorothea_hs.mat'];
        fname=fullfile(pw1,'..','resources','DoRothEA_TF_Target_DB','dorothea_hs.mat');
    case {'mm','mouse'}
        %fname=[wrkpth 'dorothea_mm.mat'];
        fname=fullfile(pw1,'..','resources','DoRothEA_TF_Target_DB','dorothea_mm.mat');
    otherwise
        error('TF database is not available for the species.');
end
fprintf('\nReading ... %s.\n',fname);
load(fname,'T');
T=T(T.mor>0,:);

if typeid==0
    score=[];
    return;
end
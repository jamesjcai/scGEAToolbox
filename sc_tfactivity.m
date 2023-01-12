function [s,tflist,gcommon]=sc_tfactivity(X,g,T,species)
% The activity level of a transcription factor (TF) in a given cell is the
% extent to which it is exerting its regulatory potential on its target 
% genes.
% https://academic.oup.com/bioinformatics/article/37/9/1234/5949002

if nargin<4, species='hs'; end
if nargin<3
    folder=fileparts(mfilename('fullpath'));
    wrkpth=fullfile(folder,'resources',filesep,'DoRothEA_TF_Target_DB',filesep);
    switch species
        case 'hs'
            fname=[wrkpth 'dorothea_hs.mat'];
        case 'mm'
            fname=[wrkpth 'dorothea_mm.mat'];
        otherwise
            fname=[wrkpth 'dorothea_hs.mat'];
    end
    load(fname,'T');
    T=T(T.mor>0,:);    % only consider positive regulation
end

[X]=sc_norm(X);

[t]=crosstab(T.tf,T.target);   % TF-by-target regulagory relationship matrix


[gcommon,k,l]=intersect(g,T.target);
[tid,tflist]=grp2idx(T.tf);

t=t(:,tid(l));     % tf-by-gene
X=X(k,:);          % gene-by-cell


%[gid,gnlist]=grp2idx(g);

s=t*X;



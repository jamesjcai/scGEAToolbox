function [cs,tflist,gcommon]=sc_tfactivity(X,g,T,species)
% The activity level of a transcription factor (TF) in a given cell is the
% extent to which it is exerting its regulatory potential on its target 
% genes.
% https://academic.oup.com/bioinformatics/article/37/9/1234/5949002
%
% [cs,tflist]=sc_tfactivity(X,g);
% CS - is an m-by-n matrix of activity scores for m TFs and n cells.
% TFlist - list of TF genes

if nargin<2, error('USAGE: [cs,tflist]=sc_tfactivity(X,g);'); end
if nargin<4, species='hs'; end
if nargin<3 || isempty(T)
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
    [gid,gnlist]=grp2idx(T.target);
    [tid,tflist]=grp2idx(T.tf);
    t=zeros(max(tid),max(gid));
    for k=1:length(gid)
        t(tid(k),gid(k))=T.mor(k);
    end
    % T=T(T.mor>0,:);    % only consider positive regulation
    %[t]=crosstab(T.tf,T.target);   % TF-by-target regulagory relationship
    %matrix if only positive regulation
end

[X]=sc_norm(X);
[X]=log(X+1);
[~,k,l]=intersect(upper(g),upper(gnlist));


t=t(:,tid(l));     % tf-by-gene
X=X(k,:);          % gene-by-cell
cs=t*X;

if nargout>2, gcommon=g(k); end

end
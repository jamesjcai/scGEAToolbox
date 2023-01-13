function [cs,tflist,gcommon]=sc_tfactivity(X,g,T,species,methodid)
% The activity level of a transcription factor (TF) in a given cell is the
% extent to which it is exerting its regulatory potential on its target 
% genes.
% https://academic.oup.com/bioinformatics/article/37/9/1234/5949002
%
% [cs,tflist]=sc_tfactivity(X,g);
% CS - is an m-by-n matrix of activity scores for m TFs and n cells.
% TFlist - list of TF genes

if nargin<5 || isempty(methodid), methodid=1; end
if nargin<4 || isempty(species), species='hs'; end
if nargin<3 || isempty(T)
if nargin<2, error('USAGE: [cs,tflist]=sc_tfactivity(X,g);'); end    
    %folder=fileparts(mfilename('fullpath'));
    %wrkpth=fullfile(folder,'resources',filesep,'DoRothEA_TF_Target_DB',filesep);
 
    pw1=fileparts(mfilename('fullpath'));
    switch lower(species)
        case {'hs','human'}
            %fname=[wrkpth 'dorothea_hs.mat'];
            fname=fullfile(pw1,'resources','DoRothEA_TF_Target_DB','dorothea_hs.mat');
        case {'mm','mouse'}
            %fname=[wrkpth 'dorothea_mm.mat'];
            fname=fullfile(pw1,'resources','DoRothEA_TF_Target_DB','dorothea_mm.mat');
        otherwise
            error('TF database is not available for the species.');
    end
    load(fname,'T');
end



methodid=1;
% T=T(T.mor>0,:);
if methodid==1, T=T(T.mor>0,:); end
if methodid==2
    [X]=sc_norm(X);
    [X]=log(X+1);
end

[gid,gnlist]=grp2idx(T.target);
[tid,tflist]=grp2idx(T.tf);
t=zeros(max(tid),max(gid));
t(sub2ind([max(tid),max(gid)],tid,gid))=T.mor;

%t2=zeros(max(tid),max(gid));
%for k=1:length(gid)
%    t2(tid(k),gid(k))=T.mor(k);
%end
%isequal(t,t2)
% T=T(T.mor>0,:);    % only consider positive regulation
%[t]=crosstab(T.tf,T.target);   % TF-by-target regulagory relationship
%matrix if only positive regulation

[~,k,l]=intersect(upper(g),upper(gnlist));
t=t(:,gid(l));     % tf-by-gene
X=X(k,:);          % gene-by-cell
if nargout>2, gcommon=g(k); end

switch methodid
    case 1        % UCell method
        cs=zeros(size(t,1),size(X,2));
        R=tiedrank(-X);
        R(R>1500)=1500+1;        
        for k=1:size(t,1)
            idx1=t(k,:)>0;
            if tflist(k)=="ZNF445"
                n1=sum(idx1)
                gcommon(idx1)
            else
                n1=sum(idx1);
            end
            if n1>0
                u=sum(R(idx1,:))-(n1*(n1-1))/2;
                cs(k,:) = 1-u/(n1*1500);
            end
        end
        cs(cs<0)=0;
    case 2    % naive method
        cs=t*X;
end




end
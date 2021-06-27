function [T]=alona(X,genelist,clusterid,varargin)

% https://alona.panglaodb.se/
% https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041
% REF: PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data

if isempty(X) || isempty(genelist)
    k=1;
    T=table("Unknown",0,'VariableNames',...
        {sprintf('C%d_Cell_Type',k),sprintf('C%d_CTA_Score',k)});
    return;
end
if nargin<3 || isempty(clusterid)
    clusterid=ones(1,size(X,2));
end
if min(size(clusterid))~=1 || ~isnumeric(clusterid)
    error('CLUSTERID={vector|[]}');
end

   p = inputParser;
   addRequired(p,'X',@isnumeric);
   addRequired(p,'genelist',@isstring);
   addRequired(p,'clusterid',@isnumeric);   
   addOptional(p,'species',"human",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["human","mouse"]));
   addOptional(p,'organ',"all",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["all","heart","immunesystem","brain","pancreas"]));   
   parse(p,X,genelist,clusterid,varargin{:});
   species=p.Results.species;
   organ=p.Results.organ;   

oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
if strcmpi(organ,"all")
    pth=fullfile(pw1,'thirdparty','alona_panglaodb2021');
else
    pth=fullfile(pw1,'thirdparty','alona_panglaodb2021',sprintf('%s',organ));
end
cd(pth);
if issparse(X)
    try
        X=full(X);
    catch
        disp('Using sparse input--longer running time is expected.');
    end
end
% warning off
% X=sc_norm(X,"type","deseq");
% warning on
X=sc_norm(X);
% disp('Library-size normalization...done.')
X=log(X+1);
genelist=upper(genelist);

switch lower(species)
    case 'human'
        Tw=readtable('markerweight_hs.txt');
        T1=readtable('markerlist_hs.txt','ReadVariableNames',false,'Delimiter','\t');
        if exist('xxmarkerlist_hs_custom.txt','file')
            T2=readtable('xxmarkerlist_hs_custom.txt','ReadVariableNames',false,'Delimiter','\t');
        else
            T2=[];
        end
    case 'mouse'
        Tw=readtable('markerweight_mm.txt');
        T1=readtable('markerlist_mm.txt','ReadVariableNames',false,'Delimiter','\t');
        if exist('xxmarkerlist_mm_custom.txt','file')
            T2=readtable('xxmarkerlist_mm_custom.txt','ReadVariableNames',false,'Delimiter','\t');
        else
            T2=[];
        end
end

wvalu=Tw.Var2;
wgene=string(Tw.Var1);
if ~isempty(T2)
    Tm=[T1;T2];
else
    Tm=T1;
end

%celltypev=string(Tm.Var1);
%[celltypev,idx]=unique(celltypev);
%Tm=Tm(idx,:);

celltypev=string(Tm.Var1);
markergenev=string(Tm.Var2);
NC=max(clusterid);

S=zeros(length(celltypev),NC);

for j=1:length(celltypev)
    g=strsplit(markergenev(j),',');
    g=strtrim(g);
    if strlength(g(end))==0
        g=g(1:end-1);
    end
    g=upper(unique(g));
    %[~,idx]=ismember(g,genelist);
    Z=zeros(NC,1); ng=zeros(NC,1);
    for i=1:length(g)
        if any(g(i)==wgene) && any(g(i)==genelist)
            wi=wvalu(g(i)==wgene);
            for k=1:NC
                z=X(g(i)==genelist,clusterid==k);
                z=mean(z(:));
                Z(k)=Z(k)+z*wi;
                ng(k)=ng(k)+1;
            end
        end
    end
    for k=1:NC
        if ng(k)>0
            S(j,k)=Z(k)./nthroot(ng(k),3);
        else
            S(j,k)=0;
        end
    end
end
T=table();
% t=table(celltypev);
for k=1:NC
    [c,idx]=sort(S(:,k),'descend');
    T=[T,table(celltypev(idx),c,'VariableNames',...
        {sprintf('C%d_Cell_Type',k),sprintf('C%d_CTA_Score',k)})];
end
if size(T,1)>10
    T=T(1:10,:);
end
cd(oldpth);
end
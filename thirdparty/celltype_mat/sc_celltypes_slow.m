function [T]=sc_celltypes_slow(X,genelist,clusterid)

% https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041
% REF: PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data

oldpth=pwd;
% pw1=fileparts(which(mfilename));
pw1 = fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/celltype_mat');
cd(pth);

X=sc_norm(X,"type","deseq");
genelist=upper(genelist);

Tw=readtable('markerweight.txt');
wvalu=Tw.Var2;
wgene=string(Tw.Var1);

Tm=readtable('markerlist.txt','ReadVariableNames',false);
celltypev=string(Tm.Var1);
markergenev=string(Tm.Var2);

S=zeros(length(celltypev),max(clusterid));

T=table(celltypev);
for k=1:max(clusterid)
    Xk=X(:,clusterid==k);
    for j=1:length(celltypev)
        g=strsplit(markergenev(j),',');
        g=g(1:end-1);
        %[~,idx]=ismember(g,genelist);
        Z=0; ng=0;
        for i=1:length(g)
            if any(g(i)==wgene) && any(genelist==g(i))
            %if ismember(g(i),wgene) && ismember(g(i),genelist)
                wi=wvalu(g(i)==wgene);                
                z=median(Xk(genelist==g(i),:));
                Z=Z+z*wi;
                ng=ng+1;
            end
        end
        if ng>0
            S(j,k)=Z./nthroot(ng,3);
        else
            S(j,k)=0;
        end
    end    
end
[~,idx]=sort(sum(S,2),'descend');

T=[T,array2table(S)];
% [~,idx]=sort(sum(table2array(T(:,2:end)),2),'descend');
T=T(idx,:);
cd(oldpth);

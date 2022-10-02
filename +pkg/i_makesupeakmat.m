function [newX,newg]=i_makesupeakmat(sce,K)
if nargin<2, K=5; end

% [chrv,stav,endv]=pkg.i_bed2mat(sce.g);
%if isstring(stav), stav=str2double(stav); end
%if isstring(endv), endv=str2double(endv); end

p=cell2mat(strfind(sce.g,":"));
chrv=extractBefore(input,p);

uchrv=unique(chrv,'stable');
newX=[];
newg=[];

for k=1:length(uchrv)
    %uchrv(k)
    ik=chrv==uchrv(k);
    %istav=stav(ik);
    ig=sce.g(ik);
    ix=sce.X(ik,:);
    a=1:K:sum(ik);
    a=[a(1:end-1); a(2:end)-1];
    iX=zeros(size(a,2),size(ix,2));
    iG=strings(size(a,2),1);
    for kk=1:size(a,2)
        iX(kk,:)=sum(ix(a(1,kk):a(2,kk),:));
        iG(kk)=ig(a(1,kk)+2);
    end
    newX=[newX;iX];
    newg=[newg;iG];
end

end
function [x,y]=e_msigjaccd(targetg,species)

if nargin<2, species='hs'; end
[~,~,Col]=gui.i_selectMSigDBGeneSet(species,true);

targetg=upper(string(targetg));

setnames=fields(Col);
x=zeros(length(setnames),1);
y=zeros(length(setnames),1);

for k=1:length(setnames)
    glist=upper(string(Col.(setnames{k}).geneSymbols));
    x(k)=length(glist);
    y(k)=length(intersect(targetg,glist)) / length(union(targetg,glist));
end

figure; 
scatter(x,y,'o');
xlabel('MsigDB gene set size');
ylabel('Jaccard similarity');
title('Gene Set - MsigDB Comparision');
dt=datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1,setnames};
end


function txt = i_myupdatefcn1(~,event_obj,g)
    idx = event_obj.DataIndex;
    txt = strrep(g(idx),'_','\_');
end



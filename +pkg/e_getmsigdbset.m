function [C]=e_getmsigdbset
urllist={'https://scgeatool.github.io/data/msigdb/c2.all.v2023.2.Hs.json',...
'https://scgeatool.github.io/data/msigdb/c5.all.v2023.2.Hs.json'};

C=struct;
Col = webread(urllist{1});
setnames = fields(Col);
for k=1:length(setnames)
    if contains(setnames(k),["KEGG_","BIOCARTA_","REACTOME_"])
        %Col=rmfield(Col,setnames(k));
        C.(setnames{k})=Col.(setnames{k});
    end
end

Col = webread(urllist{2});
setnames = fields(Col);
for k=1:length(setnames)
    if contains(setnames(k),["GOBP_","GOMF_"])
        %Col=rmfield(Col,setnames(k));
        C.(setnames{k})=Col.(setnames{k});
    end
end

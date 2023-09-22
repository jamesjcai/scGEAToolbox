function [GeneSetLibrary] = i_get_enrichr_libraries

value = webread('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=meta');
% value = jsondecode(txt);
n = length(value.libraries);
GeneSetLibrary = strings(n, 1);
for k = 1:n
    GeneSetLibrary(k) = string(value.libraries{k}.libraryName);
end
GeneSetLibrary = sort(GeneSetLibrary);


% ref: https://github.com/zqfang/GSEApy/blob/149f3ec0fabaae1ceefc1eec4a73ec423fd41cb1/gseapy/parser.py#L150

% gui.i_selmultidlg(pkg.i_get_enrichr_libraries)
end
function [GeneSetLibrary] = i_get_enrichr_libraries

value = webread('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=meta');
% value = jsondecode(txt);
n = length(value.libraries);
GeneSetLibrary = strings(n, 1);
for k = 1:n
    try
        GeneSetLibrary(k) = string(value.libraries{k}.libraryName);
    catch ME
        % warning(ME.message);
    end
end
GeneSetLibrary(GeneSetLibrary == "") = [];
GeneSetLibrary = sort(GeneSetLibrary);


% ref: https://github.com/zqfang/GSEApy/blob/149f3ec0fabaae1ceefc1eec4a73ec423fd41cb1/gseapy/parser.py#L150

% gui.i_selmultidialog(pkg.i_get_enrichr_libraries)
end
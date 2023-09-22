function [txt] = i_getsctypemarkers
txt = '';
mfolder = fileparts(mfilename('fullpath'));
infile = fullfile(mfolder, '..', 'resources', 'ScTypeDB_full.xlsx');
if exist(infile, "file")
    T = readtable(infile);
end
utissuelist = unique(T.tissueType);
[indx1, tf1] = listdlg('PromptString', ...
    {'Select Tissue Type(s):'}, ...
    'SelectionMode', 'single', ...
    'ListString', utissuelist, 'ListSize', [220, 300]);

if tf1 ~= 1, return; end
selectedtissue = utissuelist(indx1);
Tm = T(ismember(T.tissueType, selectedtissue), :);
txt = '';
for k = 1:size(Tm, 1)
    txt = sprintf('%s\n%s\t%s', txt, ...
        strtrim(Tm.cellName{k}), ...
        strtrim(Tm.geneSymbolmore1{k}));
end
txt = strtrim(txt);
end
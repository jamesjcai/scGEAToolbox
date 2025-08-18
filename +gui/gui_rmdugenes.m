function [needupdate, sce] = gui_rmdugenes(sce, parentfig)

    if nargin<2, parentfig = []; end    
    needupdate = false;
    hasDuplicates = numel(unique(sce.g)) < numel(sce.g);
    if ~hasDuplicates, return; end

    messagetext = sprintf("Method 1. Appending an underscore and a " + ...
        "number to names of duplicates\nMethod 2. Keeping " + ...
        "only the first occurrence of duplicates\nMethod " + ...
        "3. Collapsing/merging duplicates by summing their " + ...
        "expression rows\nCancel - Take no action");
    answer = gui.myQuestdlg(parentfig, messagetext, ...
        "Duplicate Genes Detected",{'Method 1', 'Method 2', ...
        'Method 3'}, {'Method 1'}, 'warning');
    switch answer
        case 'Cancel'
            return;
        case 'Method 1'
            [sce.X, sce.g] = sc_rmdugenes(sce.X, sce.g, 1);
            needupdate = true;
        case 'Method 2'
            [sce.X, sce.g] = sc_rmdugenes(sce.X, sce.g, 2);
            needupdate = true;
        case 'Method 3'
            [sce.X, sce.g] = sc_rmdugenes(sce.X, sce.g, 3);
            needupdate = true;            
        otherwise
            return;
    end
end

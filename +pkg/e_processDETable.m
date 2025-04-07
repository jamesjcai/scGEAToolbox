function [Tup, Tdn, paramset] = e_processDETable(T, paramset, parentfig)

if nargin<3, parentfig = []; end
Tup=[];
Tdn=[];
if nargin < 2, paramset = []; end
if isempty(paramset)
    [paramset] = gui.i_degparamset(false, parentfig);
end
if isempty(paramset)
   return;
else
    mindiffpct = paramset{1};
    minabsolfc = paramset{2};
    apvaluecut = paramset{3};
    sortbywhat = paramset{4};
end

a = T.(T.Properties.VariableNames{8}) - T.(T.Properties.VariableNames{7});

isok = abs(a) >= mindiffpct & abs(T.avg_log2FC) >= minabsolfc & ...
       T.p_val_adj <= apvaluecut;

fprintf(['\nDE genes with > %.2f%% difference in expression percentages, ' ...
    'abs(log2FC) >= %.2f, and adjusted P-value < %.3f are retained.\n'], ...
    mindiffpct*100,...
    minabsolfc, apvaluecut);

Tup = T(T.avg_log2FC > 0 & isok, :);
Tdn = T(T.avg_log2FC < 0 & isok, :);

if ~isempty(sortbywhat) 
    answer = sortbywhat;
else
    answer = gui.myQuestdlg(parentfig, 'Sort DE genes by adjusted P-value or fold change?','',...
        {'Adjusted P-value','Fold Change'},'Adjusted P-value');
end

switch answer
    case 'Adjusted P-value'
        Tup = sortrows(Tup, 'abs_log2FC', 'descend');
        Tup = sortrows(Tup, 'p_val_adj', 'ascend');        
        Tdn = sortrows(Tdn, 'abs_log2FC', 'descend');
        Tdn = sortrows(Tdn, 'p_val_adj', 'ascend');
        disp('DE genes are sorted by adjusted P-value.');
    case 'Fold Change'
        Tup = sortrows(Tup, 'p_val_adj', 'ascend');
        Tup = sortrows(Tup, 'abs_log2FC', 'descend');
        Tdn = sortrows(Tdn, 'p_val_adj', 'ascend');
        Tdn = sortrows(Tdn, 'abs_log2FC', 'descend');
        disp('DE genes are sorted by absolute fold change (FC).');
    otherwise
        gui.myHelpdlg(parentfig, 'Keep DE gene tables unsorted.');
end

function [Tup, Tdn] = e_processDETable(T, isatac, mindiffpct, ...
                            minabsolfc, apvaluecut)

Tup=[];
Tdn=[];
if nargin < 2, isatac = false; end
if nargin < 5
    mindiffpct = [];
    minabsolfc = [];
    apvaluecut = [];
end

if isempty(mindiffpct) || isempty(minabsolfc) || isempty(apvaluecut)
    definput = {'0.05', '1.0', '0.01'};
    prompt = {'Min. abs(diff(pct)):', ...
              'Min. abs(log2(FC)):', ...
              'Adjusted P-value cutoff:'};
    dlgtitle = 'DE Result Filter';
    dims = [1, 80];
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer), return; end
    try
        mindiffpct = str2double(answer{1});
        minabsolfc = str2double(answer{2});
        apvaluecut = str2double(answer{3});
        assert((mindiffpct > 0) && (mindiffpct < 1));
        assert((minabsolfc >= 1) && (minabsolfc < 100));
        assert((apvaluecut >= 0.0) && (apvaluecut <= 1.0));
    catch
        errordlg('Invalid input.');
        return;
    end
end

%a=T.pct_2-T.pct_1;
a = T.(T.Properties.VariableNames{8}) - T.(T.Properties.VariableNames{7});
if isatac
    isok = T.p_val_adj <= 0.1;
else
    isok = (abs(a) >= mindiffpct | abs(T.avg_log2FC) >= minabsolfc) & ...
               T.p_val_adj <= apvaluecut;
end

disp('DE genes with > 5% difference in expression percentages, abs(log2FC) >= 1, and adjusted P-value < 0.01 are retained.');

Tup = T(T.avg_log2FC > 0 & isok, :);
Tdn = T(T.avg_log2FC < 0 & isok, :);

answer = questdlg('Sort DE genes by adjusted P-value or fold change?','',...
    'Adjusted P-value','Fold Change','Adjusted P-value');

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
        waitfor(helpdlg('Keep DE gene tables unsorted.'));
end
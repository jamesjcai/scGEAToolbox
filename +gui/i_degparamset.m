function [paramset] = i_degparamset

paramset = [];
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
    answer = questdlg('Sort DE genes by adjusted P-value or fold change?','',...
        'Adjusted P-value','Fold Change','Adjusted P-value');
switch answer
    case 'Adjusted P-value'
        sortbywhat = 'Adjusted P-value';
    case 'Fold Change'
        sortbywhat = 'Fold Change';
        % disp('DE genes are sorted by absolute fold change (FC).');
    otherwise
        sortbywhat = [];
end    
paramset = {mindiffpct, minabsolfc, apvaluecut, sortbywhat};

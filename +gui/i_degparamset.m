function [paramset] = i_degparamset(nogui, parentfig)

if nargin<2, parentfig = []; end

if nargin<1, nogui=false; end
% preftagname ='scimilmodelpath'
% preftagname ='openscedlgindex';
preftagname ='degtestparamset';
defaultset = getpref('scgeatoolbox', preftagname, {0.05, 1.0, 0.01, 'Adjusted P-value'});

if nogui
    paramset = {0.05, 1.0, 0.01, 'Adjusted P-value'};   
else
    paramset = [];
    definput = {num2str(defaultset{1}), num2str(defaultset{2}), num2str(defaultset{3})};
    % definput = {'0.05', '1.0', '0.01'};
    prompt = {'Min. abs(diff(pct)): e.g., 0.05=5% (default)', ...
              'Min. abs(log2(FC)): e.g., 1.0=2x (default), 0.59=1.5x, 0.26=1.2x, 1.5=2.83x:', ...
              'Adjusted P-value cutoff: e.g., 0.01 (default)'};
    dlgtitle = 'DE Result Filter';
    dims = [1, 80];
    
    if gui.i_isuifig(parentfig)
        answer = gui.myInputdlg(prompt, dlgtitle, definput, parentfig);
    else
        answer = inputdlg(prompt, dlgtitle, dims, definput);
    end
    
    if isempty(answer), return; end
    try
        mindiffpct = str2double(answer{1});
        minabsolfc = str2double(answer{2});
        apvaluecut = str2double(answer{3});
        assert((mindiffpct >= 0) && (mindiffpct <= 1));
        assert((minabsolfc >= 0) && (minabsolfc <= 100));
        assert((apvaluecut >= 0) && (apvaluecut <= 1));
    catch
        gui.myErrordlg(parentfig, 'Invalid input.');
        return;
    end
    answer = gui.myQuestdlg(parentfig, ...
        'Sort DE genes by adjusted P-value or fold change?','',...
        {'Adjusted P-value','Fold Change'}, defaultset{4});
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
    setpref('scgeatoolbox', preftagname, paramset);
end

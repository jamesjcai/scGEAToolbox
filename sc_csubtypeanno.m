function [sce] = sc_csubtypeanno(sce, cell_type_target, speciestag)
    if nargin < 3 || isempty(speciestag), speciestag='human'; end


    pw1 = fileparts(mfilename('fullpath'));
    pth1 = fullfile(pw1, '+run', 'thirdparty', 'alona_panglaodb','marker_hs.mat');
    load(pth,'Tm');
    if ~ismember(upper(cell_type_target), upper(string(Tm.Var1)))
        error('Target cell type is not an available primary cell type.');
    end
    idx = upper(string(Tm.Var1)) == upper(cell_type_target);
    primarymarkers = Tm.Var2{idx};


    pth2 = fullfile(pw1, 'resources', 'cellsubtypes.xlsx');
    T = readtable(pth2);
    if ~ismember(upper(cell_type_target), upper(string(unique(T.CellType)))
        error('Target cell type is not a primary cell type in cellsubtype.xlsx.');
    end    
    idx = upper(string(T.CellType)) == upper(cell_type_target);
    if ~any(idx), error('SC_CSUBTYPEANNO: Runtime Error.'); end
    T = T(idx,:);

    Tm = T(:,2:3);

    

    s = upper(string(Tm.PositiveMarkers));
    S = [];
    for k = 1:length(s)
        sk = string(primarymarkers) +","+s(k);
        a = strsplit(sk, ',');
        a = strtrim(a);
        if strlength(a(end)) == 0 || isempty(a(end))
            a = a(1:end-1);
        end
        
        S = [S, a];
    end
    
    %%
    N = length(S);
    t = tabulate(S);
    f = cell2mat(t(:, 3));
    if max(f) - min(f) < eps
        w = ones(N, 1);
    else
        w = 1 + sqrt((max(f) - f)/(max(f) - min(f)));
    end
    genelist = string(t(:, 1));
    Tw = table(genelist, w);
    Tw.Properties.VariableNames = {'Var1', 'Var2'};
    


end




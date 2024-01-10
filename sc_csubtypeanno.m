function [sce] = sc_csubtypeanno(sce, cell_type_target, speciestag)

    if nargin < 3 || isempty(speciestag), speciestag='human'; end


    pw1 = fileparts(mfilename('fullpath'));
    [pmarkerstr] = in_getprimarymarkers(pw1, cell_type_target);

    pth2 = fullfile(pw1, 'resources', 'cellsubtypes.xlsx');
    T = readtable(pth2);
    if ~ismember(upper(cell_type_target), upper(string(unique(T.CellType))))
        error('Target cell type is not a primary cell type in cellsubtype.xlsx.');
    end    
    idx = upper(string(T.CellType)) == upper(cell_type_target);
    if ~any(idx), error('SC_CSUBTYPEANNO: Runtime Error.'); end
    T = T(idx,:);

    Tm = T(:,2:3);
    [Tm] = in_addprimarymarkers(Tm, pmarkerstr);

    s = upper(string(Tm.PositiveMarkers));
    S = [];
    for k = 1:length(s)
        sk = s(k);
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
    



            wvalu = Tw.Var2;
            wgene = string(Tw.Var1);
            celltypev = string(Tm.Var1);
            markergenev = string(Tm.Var2);


            sce2 = sce.selectcells(sce.c_cell_type_tx == cell_type_target);

            [Tct] = pkg.e_determinecelltype(sce2, ptsSelected, wvalu, ...
                wgene, celltypev, markergenev);            

            ctxt = Tct.C1_Cell_Type;
            if length(ctxt) > 1
                ctxt = Tct.C1_Cell_Type{1};
            end
            

end



    function [primarymarkerstr] = in_getprimarymarkers(pw1, cell_type_target)
        pth1 = fullfile(pw1, '+run', 'thirdparty', 'alona_panglaodb','marker_hs.mat');
        load(pth1,'Tm');
        if ~ismember(upper(cell_type_target), upper(string(Tm.Var1)))
            error('Target cell type is not an available primary cell type.');
        end
        idx = upper(string(Tm.Var1)) == upper(cell_type_target);
        primarymarkerstr = Tm.Var2{idx};
        primarymarkerstr = strtrim(primarymarkerstr);
        primarymarkerstr = erase(primarymarkerstr," ");
        primarymarkerstr = strip(primarymarkerstr,'right',',');
    end

    function [Tm] = in_addprimarymarkers(Tm, pmarkerstr)
        for k=1:size(T,1)
            a = string(Tm.PositiveMarkers{k});
            a = strtrim(a);
            a = erase(a," ");
            a = strip(a,'right',',');
            Tm.PositiveMarkers{k} = char(string(a)+","+pmarkerstr);
        end
    end

function [sce] = sc_csubtypeanno(sce, cell_type_target, formatid, speciestag)
    if nargin < 4 || isempty(speciestag)
        speciestag = 'human';
    end

    if nargin < 3 || isempty(formatid)
        formatid = 0;
    end

    selectedidx = sce.c_cell_type_tx == cell_type_target;

    if ~any(selectedidx)
        error('SCE.C_CELLTYPE_TXT does not contain the target cell type.');
    end

    pw1 = fileparts(mfilename('fullpath'));
    pth2 = fullfile(pw1, 'assets', 'PanglaoDB', 'cellsubtypes.xlsx');

    T = readtable(pth2);

    if ~ismember(upper(cell_type_target), upper(string(unique(T.CellType))))
        error('Target cell type is not a primary cell type in cellsubtype.xlsx.');
    end

    idx = upper(string(T.CellType)) == upper(cell_type_target);

    if ~any(idx)
        error('SC_CSUBTYPEANNO: Runtime Error.');
    end

    T = T(idx, :);

    [pmarkerstr] = in_getprimarymarkers(pw1, cell_type_target);
    Tm = T(:, 2:3);
    [Tm] = in_addprimarymarkers(Tm, pmarkerstr);

    Tw = pkg.e_markerweight(Tm);

    wvalu = Tw.Var2;
    wgene = string(Tw.Var1);
    celltypev = string(Tm.SubType);
    markergenev = string(Tm.PositiveMarkers);

    sce2 = copy(sce);
    sce2 = sce2.selectcells(selectedidx); %#OK
    sce2 = sce2.embedcells('tsne3d', true, true, 3);
    sce2 = sce2.clustercells([], [], true);

    [c, cL] = findgroups(string(sce2.c_cluster_id));

    for ik = 1:max(c)
        ptsSelected = c == ik;
        [Tct] = pkg.e_determinecelltype(sce2, ptsSelected, wvalu, ...
                wgene, celltypev, markergenev);

        ctxt = Tct.C1_Cell_Type{1};
        cL{ik} = ctxt;
    end

    sce2.c_cell_type_tx = string(cL(c));

    sce.c_cell_type_tx(selectedidx) = in_formatsubtype(sce2.c_cell_type_tx, ...
        cell_type_target, formatid);
end

function [a] = in_formatsubtype(a, b, formatid)
    switch formatid
        case 0
            return;
        case 1
            for k = 1:length(a)
                a(k) = sprintf('%s_{%s}', b, a(k));
            end
        case 2
            for k = 1:length(a)
                a(k) = sprintf('%s (%s)', b, a(k));
            end
    end
end

function [primarymarkerstr] = in_getprimarymarkers(pw1, cell_type_target)
    pth1 = fullfile(pw1, 'external', 'fun_alona_panglaodb', 'marker_hs.mat');
    load(pth1, 'Tm');

    if ~ismember(upper(cell_type_target), upper(string(Tm.Var1)))
        error('Target cell type is not an available primary cell type.');
    end

    idx = upper(string(Tm.Var1)) == upper(cell_type_target);
    primarymarkerstr = Tm.Var2{idx};
    primarymarkerstr = strtrim(primarymarkerstr);
    primarymarkerstr = erase(primarymarkerstr, " ");
    primarymarkerstr = strip(primarymarkerstr, 'right', ',');
end

function [Tm] = in_addprimarymarkers(Tm, pmarkerstr)
    for k = 1:size(Tm, 1)
        a = string(Tm.PositiveMarkers{k});
        a = strtrim(a);
        a = erase(a, " ");
        a = strip(a, 'right', ',');
        Tm.PositiveMarkers{k} = char(string(a) + "," + pmarkerstr);
    end
end

%{
function [sce] = sc_csubtypeanno(sce, cell_type_target, formatid, speciestag)

    if nargin < 4 || isempty(speciestag), speciestag='human'; end
    if nargin < 3 || isempty(formatid), formatid = 0; end

    selectedidx = sce.c_cell_type_tx == cell_type_target;
    if ~any(selectedidx)
        error('SCE.C_CELLTYPE_TXT does not contain the target cell type.');
    end


    pw1 = fileparts(mfilename('fullpath'));
    

    pth2 = fullfile(pw1, 'assets', 'PanglaoDB', 'cellsubtypes.xlsx');
    T = readtable(pth2);
    if ~ismember(upper(cell_type_target), upper(string(unique(T.CellType))))
        error('Target cell type is not a primary cell type in cellsubtype.xlsx.');
    end

    idx = upper(string(T.CellType)) == upper(cell_type_target);
    if ~any(idx), error('SC_CSUBTYPEANNO: Runtime Error.'); end
    T = T(idx,:);


    [pmarkerstr] = in_getprimarymarkers(pw1, cell_type_target);
    Tm = T(:,2:3);
    [Tm] = in_addprimarymarkers(Tm, pmarkerstr);

    [Tw] = pkg.e_markerweight(Tm);
    % s = upper(string(Tm.PositiveMarkers));
    % S = [];
    % for k = 1:length(s)
    %     sk = s(k);
    %     a = strsplit(sk, ',');
    %     a = strtrim(a);
    %     if strlength(a(end)) == 0 || isempty(a(end))
    %         a = a(1:end-1);
    %     end        
    %     S = [S, a];
    % end
    % 
    % %%
    % N = length(S);
    % t = tabulate(S);
    % f = cell2mat(t(:, 3));
    % if max(f) - min(f) < eps
    %     w = ones(N, 1);
    % else
    %     w = 1 + sqrt((max(f) - f)/(max(f) - min(f)));
    % end
    % genelist = string(t(:, 1));
    % Tw = table(genelist, w);
    % Tw.Properties.VariableNames = {'Var1', 'Var2'};
    
    wvalu = Tw.Var2;
    wgene = string(Tw.Var1);
    celltypev = string(Tm.SubType);
    markergenev = string(Tm.PositiveMarkers);


    sce2 = sce.selectcells(selectedidx); %#OK
    sce2 = sce2.embedcells('tsne3d', true, true, 3);
    sce2 = sce2.clustercells([], [], true);


    [c, cL] = findgroups(string(sce2.c_cluster_id));
    for ik = 1:max(c)
        ptsSelected = c == ik;
        [Tct] = pkg.e_determinecelltype(sce2, ptsSelected, wvalu, ...
                wgene, celltypev, markergenev);
        ctxt = Tct.C1_Cell_Type{1};
        cL{ik} = ctxt;
    end
    sce2.c_cell_type_tx = string(cL(c));

    sce.c_cell_type_tx(selectedidx) = in_formatsubtype(sce2.c_cell_type_tx, ...
        cell_type_target, formatid);

end


    function [a]=in_formatsubtype(a, b, formatid)
        switch formatid
            case 0
                return;
            case 1
                for k = 1:length(a)
                    a(k) = sprintf('%s_{%s}',b,a(k));
                end
            case 2
                for k = 1:length(a)
                    a(k) = sprintf('%s (%s)',b,a(k));
                end                
        end
    end


%}


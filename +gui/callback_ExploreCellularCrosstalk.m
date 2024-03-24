function callback_ExploreCellularCrosstalk(src, ~)

FigureHandle = src.Parent.Parent;
if ~gui.gui_showrefinfo('talklr [DOI:10.1101/2020.02.01.930602]'), return; end

    answer = questdlg('This function is based on an unpublished method [DOI:10.1101/2020.02.01.930602]. Continue?');
        if ~strcmp(answer, 'Yes'), return; end
        
        sce = guidata(FigureHandle);
        if isempty(sce.c_cell_type_tx) || numel(unique(sce.c_cell_type_tx)) < 2
            if ~isempty(sce.c_cluster_id) && numel(unique(sce.c_cluster_id)) > 1
                answer = questdlg(sprintf('Cell type (C_CELL_TYPE_TX) is undefined.\nWould you like to use cluster id (C_CLUSTER_ID) to define cell groups?'));
                switch answer
                    case 'Yes'
                        sce.c_cell_type_tx = strcat('Goup', string(sce.c_cluster_id));
                    otherwise
                        return;
                end
            else
                warndlg('Cell type is undefined (SCE.C_CELL_TYPE_TX is empty)');
                return;
            end
        end

        [c, cL] = grp2idx(sce.c_cell_type_tx);
        [idx] = gui.i_selmultidlg(cL, [], FigureHandle);
        if isempty(idx), return; end
        if numel(idx) < 2
            warndlg('Need at least 2 cell types');
            return;
        end
        selected = ismember(c, idx);


        pw = fileparts(mfilename('fullpath'));
        dbfile = fullfile(pw, '..', 'resources', 'Ligand_Receptor.mat');
        load(dbfile, 'ligand', 'receptor', 'T');


        fw = gui.gui_waitbar;
        sce = sce.selectcells(selected);
        [OUT, ~] = run.mt_talklr(sce);
        gui.gui_waitbar(fw);

        n = length(OUT.ligandok);
        if n == 0
            warndlg('Not detected.');
        end

        labels = {'Save OUT to variable named:'};
        vars = {'OUT'};
        values = {OUT};

        if ~(ismcc || isdeployed)
            [f, ft] = export2wsdlg(labels, vars, values);
            waitfor(f);
            if ft
                disp('Run gui.i_crosstalkgraph(OUT,k) to plot crosstalk graph for ligand-receptor pair k.')
            end
        end

            listitems = cell(n, 1);
            for k = 1:n
                listitems{k} = sprintf('%s -> %s (KL = %.2f)', ...
                    OUT.ligandok(k), OUT.receptorok(k), ...
                    OUT.KL(k));
            end
            gui.i_exporttable(T, true, 'Tccrosstalk', 'CrossTalkTable');            

            i_displyres(listitems);


            function i_displyres(listitems)
                [indx2, tf2] = listdlg('PromptString', ...
                    {'Select ligand-receptor pairs to plot'}, ...
                    'SelectionMode', 'single', 'ListString', listitems, ...
                    'ListSize', [210, 300]);
                if tf2 == 1
                    kk = indx2;
                    [y1, idx1] = ismember(upper(OUT.ligandok(kk)), upper(sce.g));
                    [y2, idx2] = ismember(upper(OUT.receptorok(kk)), upper(sce.g));
                    if y1 && y2
                        hFig = figure;
                        hFig.Position(3) = hFig.Position(3) * 2.2;

                        %subplot(1, 2, 1)
                        nexttile
                        sc_scattermarker(sce.X, sce.g, sce.s, sce.g(idx1), 1, [], false);
                        title(sce.g(idx1));

                        %subplot(1, 2, 2)
                        nexttile
                        sc_scattermarker(sce.X, sce.g, sce.s, sce.g(idx2), 1, [], false);
                        title(sce.g(idx2));
                        
                        tb = uitoolbar(hFig);
                        pt5pickcolr = uipushtool(tb, 'Separator', 'off');
                        [img, map] = imread(fullfile(fileparts(mfilename('fullpath')), ...
                            '../', 'resources', 'fvtool_fdalinkbutton.gif')); % plotpicker-pie
                        % map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;  % Convert white pixels => transparent background
                        ptImage = ind2rgb(img, map);
                        pt5pickcolr.CData = ptImage;
                        pt5pickcolr.Tooltip = 'Link subplots';
                        pt5pickcolr.ClickedCallback = @gui.i_linksubplots;

                        nexttile
                        
                        ligand_mat = OUT.ligand_mat;
                        receptor_mat = OUT.receptor_mat;
                        ligandok = OUT.ligandok;
                        receptorok = OUT.receptorok;
                    
                        a = ligand_mat(k, :);
                        b = receptor_mat(k, :);
                        m = (a' * b) .* ((a > 0)' * (b > 0));
                        m = m ./ sum(m(m > 0));
                    
                        % figname = sprintf('%s (ligand) -> %s (receptor)', ...
                        %     ligandok(k), receptorok(k));
                        G = pkg.i_makegraph(m, OUT.cL);
                        p = plot(G);
                        cc = repmat([0, 0.4470, 0.7410], G.numedges, 1);
                        cc(G.Edges.Weight < 0, :) = repmat([0.8500, 0.3250, 0.0980], ...
                            sum(G.Edges.Weight < 0), 1);
                        p.EdgeColor = cc;
                        if ~isempty(G.Edges.Weight)
                            w = 3;
                            G.Edges.LWidths = abs(w*G.Edges.Weight/max(G.Edges.Weight));
                            p.LineWidth = G.Edges.LWidths;
                        end    
                        % sc_grnview(m, OUT.cL, figname);

                    end

                    % gui.i_crosstalkgraph(OUT, kk);
                    i_displyres(listitems);
                elseif tf2 == 0
                    return;
                end
            end


        end

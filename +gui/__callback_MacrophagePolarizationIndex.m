        if indx == 5
            hc = colorbar;
            hc.Label.String = ttxt;
        else
            colorbar off;
        end
        [c, cL] = grp2idx(ci);
        sce.c = ci;
        guidata(FigureHandle, sce);
        
        if indx == 5
            answer = questdlg('Compare G1/S/G2M ratios between classes?');
            if ~isequal(answer, 'Yes')
                return
            end
            
            listitems = {'Cluster ID', 'Batch ID', ...
                'Cell Type'};
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select statistics', '', ''}, ...
                'SelectionMode', 'single', ...
                'ListString', listitems);
            if tf2 == 1
                switch indx2
                    case 1 % cluster id
                        thisc = sce.c_cluster_id;
                    case 2 % batch id
                        thisc = sce.c_batch_id;
                    case 3 % cell type
                        thisc = sce.c_cell_type_tx;
                end
            else
                return
            end
            if isempty(thisc)
                errordlg("Undefined classification");
                return
            end
            try
                [A, ~, ~, l] = crosstab(sce.c_cell_cycle_tx, thisc);
                B = A ./ sum(A);
                figure;
                bar(B', 'stacked');
                ylabel(sprintf('%s|', string(l(1:3, 1))));
            catch ME
                rethrow(ME);
            end
        end
        
        
% uimenu(m,'Text','Macrophage Polarization Index...',...
%     'Callback',@callback_MacrophagePolarizationIndex);


function callback_MacrophagePolarizationIndex(src,~)

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    helpdlg("Function is under development.");
    return;

    scorename="Macrophage_Polarization_Index";
    cs=pkg.e_cellscores(sce.X,sce.g,scorename);
    figure;
    gui.i_stemscatter(sce.s,cs);
    zlabel('Score Value')
    title(scorename)
    
    labels = {'Save score values to variable named:'}; 
    vars = {scorename};
    values = {cs};
    export2wsdlg(labels,vars,values);
end
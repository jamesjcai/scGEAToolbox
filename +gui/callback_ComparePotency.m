function callback_ComparePotency(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);   
    if isempty(sce.list_cell_attributes)
        warndlg('No data sce.list_cell_attributes')
        return;
    end
    if ~isempty(sce.list_cell_attributes{2}) && ~isempty(sce.c_cell_type_tx)
        x=sce.list_cell_attributes{2};
        y=sce.c_cell_type_tx;
    figure;
    pkg.i_violinplot_groupordered(x,y);
    if ~isempty(sce.c_batch_id)
        figure;
        h1=subplot(1,2,1);
        i=sce.c_batch_id==1;
        pkg.i_violinplot_groupordered(x(i),y(i));
        [aabb]=ylim();
        h2=subplot(1,2,2);
        i=sce.c_batch_id==2;
        pkg.i_violinplot_groupordered(x(i),y(i));
        [ccdd]=ylim();
        ylim(h1,[min(aabb(1),ccdd(1)) max(aabb(2),ccdd(2))]);
        ylim(h2,[min(aabb(1),ccdd(1)) max(aabb(2),ccdd(2))]);
    end
    end
end

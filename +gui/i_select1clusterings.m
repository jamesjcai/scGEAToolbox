function [thisc,clable]=i_select1clusterings(sce)
thisc=[];
clable='';

listitems={''};
methodtagv=fieldnames(sce.struct_cell_clusterings);
for k=1:length(methodtagv)
    methodtag=methodtagv{k};
    if ~isempty(sce.struct_cell_clusterings.(methodtag))
        listitems=[listitems,methodtag];
    end
end
listitems(1)=[];
if isempty(listitems), return; end

[indx2,tf2] = listdlg('PromptString',...
    {'Select clustering variable:'},...
     'SelectionMode','single','ListString',listitems);
if tf2==1
    clable=listitems{indx2};
    thisc=sce.struct_cell_clusterings.(clable);
end

end
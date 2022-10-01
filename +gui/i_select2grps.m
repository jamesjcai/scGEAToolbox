function [i1,i2,cL1,cL2]=i_select2grps(sce)
i1=0; i2=0;
cL1=[]; cL2=[];
% internal function called by callback_DEGene2Groups
% listitems={'Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};
% [indx2,tf2] = listdlg('PromptString',...
%     {'Select statistics','',''},...
%      'SelectionMode','single','ListString',listitems);
% if tf2==1
%     switch indx2
%         case 1 % cluster id
%             thisc=sce.c_cluster_id;
%         case 2 % batch id
%             thisc=sce.c_batch_id;
%         case 3 % cell type
%             thisc=sce.c_cell_type_tx;                
%         case 4 % cell cycle
%             thisc=sce.c_cell_cycle_tx;
%     end
% else
%     return;
% end 


[thisc,~]=gui.i_select1class(sce);
if isempty(thisc)
    %errordlg('Undefined');
    return;
end
if numel(unique(thisc))==1
    warndlg("Cannot compare with an unique group");
    return;
end

[ci,cLi]=grp2idx(thisc);
listitems=sort(string(cLi));
n=length(listitems);
if n<2
    errordlg('Need at least two groups.');
    return;
end
[indxx,tfx] = listdlg('PromptString',{'Select two groups:'},...
    'SelectionMode','multiple',...
    'ListString',listitems,...
    'InitialValue',[n-1 n]);
if tfx==1
    if numel(indxx)~=2
        errordlg('Please select 2 groups');
        return;
    end
    [y1,idx1]=ismember(listitems(indxx(1)),cLi);
    [y2,idx2]=ismember(listitems(indxx(2)),cLi);
    assert(y1&y2);
    i1=ci==idx1;
    i2=ci==idx2;
    cL1=cLi(idx1);
    cL2=cLi(idx2);

    %i1=ismember(ci,indxx(1));
    %i2=ismember(ci,indxx(2));
    %cL1=cLi(indxx(1));
    %cL2=cLi(indxx(2));
else
    return;
end
end
function [indx1,species]=i_selgenecollection

species=[];
indx1=[];
selitems={'MSigDB Molecular Signatures',...
          'DoRothEA TF Targets Expression',...
        'Predefined Gene Collections'};
    [indx1,tf1]=listdlg('PromptString',...
        'Select a gene set collection.',...
        'SelectionMode','single','ListString',selitems, ...
        'ListSize', [220, 300]);
    if tf1~=1, return; end

    if indx1==1
        species = gui.i_selectspecies(2);
    else
        species='human';
    end
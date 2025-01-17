function [selecteditem, speciesid] = i_selgenesetcollection

speciesid = [];
selecteditem = [];

%MSigDB Molecular Signatures
%PanglaoDB Cell Type Markers
%DoRothEA TF Targets
%Custom Gene Sets

selitems={'MSigDB Molecular Signatures',...
          'DoRothEA TF Targets',...
          'PanglaoDB Cell Type Markers',...
          'Predefined Custom Gene Sets',...
          '---------------------------',...
          'Define a New Score...'};

    [selectedindx, tf1]=listdlg('PromptString',...
        'Select a gene set collection or define a new gene set.',...
        'SelectionMode','single','ListString',selitems, ...
        'ListSize', [220, 300]);
    if tf1~=1, return; end
    
    selecteditem = selitems{selectedindx};
    if selectedindx == 5, return; end
    if selectedindx == 4, speciesid = 'human'; return; end
    speciesid = gui.i_selectspecies(2);
end

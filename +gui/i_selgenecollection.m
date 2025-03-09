function [indx1, species] = i_selgenecollection

%see also: i_selectcellscore

%MSigDB Molecular Signatures
%PanglaoDB Cell Type Markers
%DoRothEA TF Targets
%Custom Gene Sets

species=[];

selitems={'MSigDB Molecular Signatures',...
          'DoRothEA TF Targets',...
          'Custome Gene Sets'};

    [indx1, tf1]=listdlg('PromptString',...
        'Select a gene set collection.',...
        'SelectionMode','single','ListString',selitems, ...
        'ListSize', [220, 300]);
    if tf1~=1, return; end

    if indx1==1
        species = gui.i_selectspecies(2, false, FigureHandle);
    else
        species='human';
    end
function [indx1, species] = i_selgenecollection(parentfig)
if nargin < 1, parentfig = []; end
%see also: i_selectcellscore

%MSigDB Molecular Signatures
%PanglaoDB Cell Type Markers
%DoRothEA TF Targets
%Custom Gene Sets

species=[];

selitems={'MSigDB Molecular Signatures',...
          'DoRothEA TF Targets',...
          'Custome Gene Sets'};

if gui.i_isuifig(parentfig)
    [indx1, tf1] = gui.myListdlg(parentfig, selitems, 'Select a gene set collection.');
else
    [indx1, tf1]=listdlg('PromptString',...
        'Select a gene set collection.',...
        'SelectionMode','single','ListString',selitems, ...
        'ListSize', [220, 300]);
end
    if tf1~=1, return; end

    if indx1==1
        species = gui.i_selectspecies(2, false, parentfig);
    else
        species='human';
    end
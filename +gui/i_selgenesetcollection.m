function [selecteditem, speciesid] = i_selgenesetcollection(parentfig)
if nargin<1, parentfig=[]; end

speciesid = [];
selecteditem = [];

% MSigDB Molecular Signatures
% PanglaoDB Cell Type Markers
% Custom Gene Sets
% Note: TF activity analysis is available via Analyze > TF Activity Analysis

selitems={'MSigDB Molecular Signatures',...
          'PanglaoDB Cell Type Markers',...
          'Predefined Custom Gene Sets',...
          '---------------------------',...
          'Define a New Score...'};

        if gui.i_isuifig(parentfig)
            [selectedindx, tf1] = gui.myListdlg(parentfig, selitems, ...
                'Select a gene set collection or define a new gene set.');
        else
            [selectedindx, tf1]=listdlg('PromptString',...
                'Select a gene set collection or define a new gene set.',...
                'SelectionMode','single','ListString',selitems, ...
                'ListSize', [220, 300]);
        end
if tf1~=1, return; end

selecteditem = selitems{selectedindx};
if selectedindx == 4, return; end              % '---------------------------' separator
if selectedindx == 3, speciesid = 'human'; return; end  % 'Predefined Custom Gene Sets'
speciesid = gui.i_selectspecies(2, false, parentfig);  % MSigDB, PanglaoDB
end

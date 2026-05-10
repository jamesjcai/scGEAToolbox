function [indx1, species] = i_selgenecollection(parentfig, preferredspecies)
if nargin < 1, parentfig = []; end
if nargin < 2, preferredspecies = []; end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end
% see also: i_selectcellscore % OK

% MSigDB Molecular Signatures
% PanglaoDB Cell Type Markers
% DoRothEA TF Targets
% Custom Gene Sets

species=[];

if isempty(preferredspecies)
    selitems = {'MSigDB Molecular Signatures (Human)', ...
                'MSigDB Molecular Signatures (Mouse)', ...
                'DoRothEA TF Targets', ...
                'Custom Gene Sets'};
else
    selitems = {'MSigDB Molecular Signatures', ...
                'DoRothEA TF Targets', ...
                'Custom Gene Sets'};
end

if gui.i_isuifig(parentfig)
    [indx1, tf1] = gui.myListdlg(parentfig, selitems, 'Select a gene set collection.');
else
    [indx1, tf1] = listdlg('PromptString', ...
        'Select a gene set collection.', ...
        'SelectionMode', 'single', 'ListString', selitems, ...
        'ListSize', [260, 300]);
end
if tf1 ~= 1, return; end

if isempty(preferredspecies)
    switch indx1
        case 1
            species = 'human';
            indx1 = 1;
        case 2
            species = 'mouse';
            indx1 = 1;
        case 3
            species = 'human';
            indx1 = 2;
        case 4
            species = 'human';
            indx1 = 3;
    end
else
    switch indx1
        case 1
            species = preferredspecies;
        otherwise
            species = 'human';
    end
end

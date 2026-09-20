function [indx1, species] = i_selgenecollection(parentfig, preferredspecies)
if nargin < 1, parentfig = []; end
if nargin < 2, preferredspecies = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
% see also: i_selectcellscore % OK

% MSigDB Molecular Signatures
% PanglaoDB Cell Type Markers
% DoRothEA TF Targets
% Custom Gene Sets

species=[];

% Each list is paired with the PKG.E_GETGENESETS option each entry maps to.
% Keeping the two side by side is the point: the lists used to be mapped by
% a switch on the list index, which happened to agree with the option
% numbers and so hid the fact that they are different things. Appending one
% entry was enough to break that - 'GlycoEnzOnto Pathways' sits at list
% index 5 in the second list, and option 5 is DoRothEA signed.
if isempty(preferredspecies)
    selitems = {'MSigDB Molecular Signatures (Human)', ...
                'MSigDB Molecular Signatures (Mouse)', ...
                'DoRothEA TF Targets', ...
                'Custom Gene Sets', ...
                'Glycobiology Gene Sets', ...
                'GlycoEnzOnto Pathways'};
    options = [1, 1, 2, 3, 4, 6];
    speciesof = {'human', 'mouse', 'human', 'human', 'human', 'human'};
else
    selitems = {'MSigDB Molecular Signatures', ...
                'DoRothEA TF Targets', ...
                'Custom Gene Sets', ...
                'Glycobiology Gene Sets', ...
                'GlycoEnzOnto Pathways'};
    options = [1, 2, 3, 4, 6];
    speciesof = {preferredspecies, 'human', 'human', 'human', 'human'};
end

if gui.i_isuifig(parentfig)
    % allowmulti = false, explicitly: myListdlg defaults it to true, and the
    % mapping below takes one index. The listdlg branch has always said
    % 'single'; this is the uifigure branch agreeing with it. Same reasoning
    % as in GUI.I_SELECT1CLASS.
    [indx1, tf1] = gui.myListdlg(parentfig, selitems, ...
        'Select a gene set collection.', [], false);
else
    [indx1, tf1] = listdlg('PromptString', ...
        'Select a gene set collection.', ...
        'SelectionMode', 'single', 'ListString', selitems, ...
        'ListSize', [260, 300]);
end
if tf1 ~= 1, return; end
if isempty(indx1) || indx1 < 1 || indx1 > numel(options)
    indx1 = [];
    return;
end

species = speciesof{indx1};
indx1 = options(indx1);

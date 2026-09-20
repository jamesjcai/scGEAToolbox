function Tnt = i_decutoffnote(Tnt, usedset)
%I_DECUTOFFNOTE Append the DE cutoffs actually applied to a Note table.
%
%   TNT = PKG.I_DECUTOFFNOTE(TNT, USEDSET) adds one row per cutoff to TNT,
%   the Note table from PKG.IN_DETABLEPROCESS. USEDSET is the fourth
%   output of PKG.E_PROCESSDETABLE. With 'Automatic' cutoffs these can
%   differ between cell types, so each workbook has to say which it used.
%   TNT is returned unchanged when USEDSET is empty or TNT has no
%   Description column.

if isempty(usedset) || ~ismember('Description', Tnt.Properties.VariableNames)
    return;
end

cutoffmode = 'Fixed';
if numel(usedset) >= 5 && ~isempty(usedset{5})
    cutoffmode = usedset{5};
end

Item = {'Cutoff mode'; 'Min. abs(diff(pct)) applied'; ...
    'Min. abs(log2FC) applied'; 'Adjusted P-value cutoff applied'};
Description = {cutoffmode; sprintf('%g', usedset{1}); ...
    sprintf('%g', usedset{2}); sprintf('%g', usedset{3})};
% The last automatic stage can add a raw P-value cutoff to break BH ties.
if numel(usedset) >= 6 && isfinite(usedset{6})
    Item{end+1} = 'Raw P-value cutoff applied';
    Description{end+1} = sprintf('%g', usedset{6});
end
Tnt = [Tnt; table(Item, Description)];
end

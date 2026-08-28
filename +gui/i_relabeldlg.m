function newvals = i_relabeldlg(parentfig, vals, dlgtitle)
% I_RELABELDLG - Edit the distinct values of a per-cell attribute.
%
%   newvals = gui.i_relabeldlg(parentfig, vals)
%   newvals = gui.i_relabeldlg(parentfig, vals, dlgtitle)
%
% Inputs:
%   parentfig - parent figure, raised again when the dialog closes (default: [])
%   vals      - per-cell values, one element per cell
%   dlgtitle  - dialog title (default: 'Edit Distinct Values')
%
% Outputs:
%   newvals - string values with the edits applied, same size and orientation as
%             vals; [] when cancelled or when nothing was changed
%
% Shows one row per distinct value with its cell count, and writes whatever is
% typed in the New Value column back to every cell carrying that value. This is
% how a long annotation is edited without ever putting one line per cell on
% screen: a 200k-cell attribute with a dozen labels is a dozen rows here, while
% the same values in a text area take minutes to load.
%
% Two distinct values can be merged by giving them the same new value.
%
% Missing values - NaN in a numeric attribute, <missing> in a string one - are
% shown as one more distinct value, under a "<missing>" placeholder, so they can
% be counted and given a real name. Leaving that row alone writes missing back
% rather than the placeholder text.
%
% see also: gui.sc_cellattribeditor, pkg.i_grp2idxsorted

if nargin < 1, parentfig = []; end
if nargin < 3 || isempty(dlgtitle), dlgtitle = 'Edit Distinct Values'; end

newvals = [];
if isempty(vals), return; end

[grp, labels] = pkg.i_grp2idxsorted(vals);
labels = string(labels(:));
grp = grp(:);

% findgroups assigns no group to a missing value, so pkg.i_grp2idxsorted hands
% back NaN for those rows. Give them a group of their own: accumarray below
% rejects NaN subscripts outright, and a NaN would index nothing when the edited
% labels are mapped back onto the cells.
ismissinggrp = isnan(grp);
missingtag = strings(0, 1);
if any(ismissinggrp)
    % Suffixed if some cell genuinely holds the text "<missing>", so the
    % placeholder row can never be confused with a real value.
    missingtag = string(matlab.lang.makeUniqueStrings("<missing>", labels));
    labels = [labels; missingtag];
    grp(ismissinggrp) = numel(labels);
end

counts = accumarray(grp, 1, [numel(labels), 1]);

dlgSize = [420, 450];
dlgPos = gui.i_centerdlgpos(parentfig, dlgSize);
d = uifigure('Name', dlgtitle, 'Position', dlgPos, ...
    'WindowStyle', 'normal', 'Visible', 'off');

uilabel(d, 'Position', [20 dlgSize(2)-40 dlgSize(1)-40 22], ...
    'Text', sprintf('%d distinct values in %d cells. Edit the New Value column:', ...
    numel(labels), numel(grp)));

ut = uitable(d, 'Position', [20 60 dlgSize(1)-40 dlgSize(2)-110]);
ut.Data = table(labels, counts, labels, ...
    'VariableNames', {'Value', 'Cells', 'NewValue'});
ut.ColumnName = {'Value', 'Cells', 'New Value'};
ut.RowName = {};
ut.ColumnWidth = {'1x', 60, '1x'};
ut.ColumnEditable = [false, false, true];

uibutton(d, 'Text', 'OK', 'Position', [60 20 80 30], ...
    'ButtonPushedFcn', @(btn, event) in_ok(d));
uibutton(d, 'Text', 'Cancel', 'Position', [160 20 80 30], ...
    'ButtonPushedFcn', @(btn, event) uiresume(d));

if ~isMATLABReleaseOlderThan('R2025a')
    try theme(d, parentfig.Theme.BaseColorStyle); catch; end
end

d.UserData = false;
drawnow;
d.Visible = 'on';

% WindowStyle='modal' is intentionally omitted, for the reason gui.myListdlg
% records: on a multi-monitor setup whose monitors differ in DPI, MATLAB
% re-centers a modal window onto the primary monitor and discards the Position
% computed above by gui.i_centerdlgpos. uiwait blocks the caller either way, so
% the dialog is functionally modal without it.
uiwait(d);

if pkg.i_isvalid(d) && d.UserData
    newlabels = strtrim(string(ut.Data.NewValue));
    if ~isequal(newlabels, labels)
        % An untouched placeholder row means the user did not name the missing
        % cells, so they stay missing instead of picking up the placeholder text.
        if ~isempty(missingtag) && newlabels(end) == missingtag
            newlabels(end) = missing;
        end
        newvals = newlabels(grp);
        newvals = reshape(newvals, size(vals));
    end
end
if pkg.i_isvalid(d)
    uiresume(d);
    delete(d);
end

if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
    figure(parentfig);
end
end

function in_ok(d)
d.UserData = true;
uiresume(d);
end

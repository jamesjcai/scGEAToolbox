function i_addharmonymenu(app)
% I_ADDHARMONYMENU - Inject the native Harmony menu item under the Edit menu.
% scgeatoolApp.mlapp is a binary App Designer file, so the menu is added at
% launch time rather than by editing the .mlapp. The function is idempotent:
% if the item already exists (e.g. later baked into the app), it does nothing.
%
% It goes under Edit rather than R Tools or Python Tools, where the two
% existing Harmony entries live, because RUN.ML_HARMONY is native MATLAB and
% needs neither installation. Its neighbours here are the other operations
% that change the data itself -- merging SCE files, renaming batch IDs --
% and merging is what creates the multi-batch data this corrects.
%
% See also GUI.CALLBACK_HARMONY, RUN.ML_HARMONY, GUI.I_ADDINFERCNVMENU.

if nargin < 1 || isempty(app), return; end
if ~isprop(app, 'EditMenu') || ~pkg.i_isvalid(app.EditMenu), return; end

menutext = 'Batch Integration (Harmony) [PMID:31740819]...';
existing = findobj(app.EditMenu, 'Type', 'uimenu', 'Text', menutext);
if ~isempty(existing), return; end

uimenu(app.EditMenu, ...
    'Text', menutext, ...
    'Separator', 'on', ...
    'Tooltip', ['Remove batch effects from the cell embedding with the ' ...
    'native MATLAB implementation of Harmony (no R or Python required)'], ...
    'MenuSelectedFcn', @(~, ~) gui.callback_Harmony(app));
end

function i_addcogapsmenu(app)
% I_ADDCOGAPSMENU - Inject the CoGAPS menu item under the R Tools menu.
% scgeatoolApp.mlapp is a binary App Designer file, so the menu is added at
% launch time rather than by editing the .mlapp. The function is idempotent:
% if the item already exists (e.g. later baked into the app), it does nothing.

if nargin < 1 || isempty(app), return; end
if ~isprop(app, 'RToolsMenu') || ~pkg.i_isvalid(app.RToolsMenu), return; end

menutext = 'NMF Pattern Discovery (CoGAPS) [PMID:37828301]...';
existing = findobj(app.RToolsMenu, 'Type', 'uimenu', 'Text', menutext);
if ~isempty(existing), return; end

uimenu(app.RToolsMenu, ...
    'Text', menutext, ...
    'Tooltip', '[PMID:37828301]', ...
    'MenuSelectedFcn', @(~, ~) gui.callback_RunCoGAPS(app));
end

function i_addtnviewmenu(app)
% I_ADDTNVIEWMENU - Inject the scTenifoldNet viewer item under Network menu.
% scgeatoolApp.mlapp is a binary App Designer file, so the menu is added at
% launch time rather than by editing the .mlapp. The function is idempotent:
% if the item already exists (e.g. later baked into the app), it does nothing.

if nargin < 1 || isempty(app), return; end
if ~isprop(app, 'NetworkMenu') || ~pkg.i_isvalid(app.NetworkMenu), return; end

menutext = 'View scTenifoldNet Result (Two GRNs + DR)...';
existing = findobj(app.NetworkMenu, 'Type', 'uimenu', 'Text', menutext);
if ~isempty(existing), return; end

uimenu(app.NetworkMenu, ...
    'Text', menutext, ...
    'Tooltip', '[PMID:33336197]', ...
    'MenuSelectedFcn', @(~, ~) gui.callback_scTenifoldNetView(app));
end

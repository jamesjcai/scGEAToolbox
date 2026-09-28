function method = i_selectgrnmethod(parentfig, prompt)
%I_SELECTGRNMETHOD Pick a gene-network method from net.grnmethods.
%   METHOD = gui.i_selectgrnmethod(parentfig) lists the methods marked
%   InMenus and returns the chosen row of net.grnmethods() as a one-row
%   table, or [] if the dialog is cancelled. Callers read METHOD.Key for
%   sc_grn, METHOD.Transform to preselect gui.i_transformx, and pass
%   METHOD to gui.i_confirmgrnrun before running it.
%
%   METHOD = gui.i_selectgrnmethod(parentfig, prompt) shows PROMPT above
%   the list.
%
%   See also net.grnmethods, gui.i_confirmgrnrun, sc_grn.

if nargin < 1, parentfig = []; end
if nargin < 2, prompt = ''; end

catalog = net.grnmethods();
catalog = catalog(catalog.InMenus, :);

method = [];
[sel, ok] = gui.myListdlg(parentfig, cellstr(catalog.Label), ...
    'Select GRN construction method', 1, false, true, [300, 450], prompt);
if ~ok || isempty(sel), return; end
method = catalog(sel, :);
end

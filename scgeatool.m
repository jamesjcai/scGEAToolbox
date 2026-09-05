function varargout = scgeatool(varargin)
% SCGEATOOL - Launch App Designer GUI
%
% Use:
%   scgeatool            % launches App Designer GUI

if ~gui.i_installed('stats'), return; end
app = scgeatoolApp(varargin{:});
try
    gui.i_addtnviewmenu(app);
catch
    % best-effort: the GUI still works without the injected viewer menu
end
% try
%     gui.i_addcogapsmenu(app);
% catch
%     % best-effort: the GUI still works without the injected CoGAPS menu
% end
% try
%     gui.i_reorganizeannotatemenu(app);
% catch
%     % best-effort: the GUI still works with the Annotate menu unreorganized
% end
% try
%     gui.i_addcompareannotmenu(app);
% catch
%     % best-effort: the GUI still works without the injected compare menu
% end
% try
%     gui.i_addsubtypemenu(app);
% catch
%     % best-effort: the GUI still works without the injected subtype menu
% end
if nargout > 0
    varargout{1} = app;
end
end
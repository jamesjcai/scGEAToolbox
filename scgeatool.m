function varargout = scgeatool(varargin)
% SCGEATOOL - Launch App Designer GUI
%
% Use:
%   scgeatool            % launches App Designer GUI

if ~gui.i_installed('stats'), return; end
app = scgeatoolApp(varargin{:});
try
    gui.i_addcogapsmenu(app);
catch
    % best-effort: the GUI still works without the injected CoGAPS menu
end
if nargout > 0
    varargout{1} = app;
end
end
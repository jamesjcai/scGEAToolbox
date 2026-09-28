function varargout = scgeatool(varargin)
% SCGEATOOL - Launch App Designer GUI
%
% Use:
%   scgeatool            % launches App Designer GUI
%
% Every menu item, label, tooltip and the whole menu arrangement live in
% scgeatoolApp.mlapp, so this launches the app and nothing else. Change a
% menu there, not by reshaping it at launch: eight places construct
% SCGEATOOLAPP directly, and a launch-time change here would miss them.

if ~gui.i_installed('stats'), return; end
app = scgeatoolApp(varargin{:});
if nargout > 0
    varargout{1} = app;
end
end

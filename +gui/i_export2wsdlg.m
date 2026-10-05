function OKPressed = i_export2wsdlg(returnfig, labels, vars, values, varargin)
% I_EXPORT2WSDLG - EXPORT2WSDLG that hands focus back to the calling window.
%
%   OKPressed = gui.i_export2wsdlg(returnfig, labels, vars, values, ...)
%
% Takes EXPORT2WSDLG's arguments after RETURNFIG (title, defaults, ...),
% waits for the dialog to close, then raises RETURNFIG. When the dialog
% closes, Windows gives focus to the MATLAB desktop instead of the window it
% was opened from, so that window ends up behind the desktop. An empty
% RETURNFIG means the window whose callback is running (GCBF).
%
% see also: export2wsdlg, gui.i_raisefig, gui.myExport2wsdlg

if isempty(returnfig)
    returnfig = gcbf;
end
% With both outputs requested, EXPORT2WSDLG blocks until the dialog closes.
[~, OKPressed] = export2wsdlg(labels, vars, values, varargin{:});
gui.i_raisefig(returnfig);
end

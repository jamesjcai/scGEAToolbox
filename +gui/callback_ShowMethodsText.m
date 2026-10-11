function callback_ShowMethodsText(src, ~)
%CALLBACK_SHOWMETHODSTEXT Show the Methods section for the analysis so far.
%
%   gui.callback_ShowMethodsText(app)
%
% File > Export Methods Text. Collects the paragraph each analysis step
% recorded on the data (SC_METHODSTEXT), prints it to the Command Window
% and shows it in a window with Copy to Clipboard and Save as Text File
% buttons (GUI.I_SHOWMETHODSTEXT). Nothing reaches the clipboard unless the
% user clicks Copy. Changes nothing.
%
% see also: sc_methodstext, gui.i_showmethodstext, pkg.i_methodsreport

[FigureHandle, sce] = gui.gui_getfigsce(src);
txt = sc_methodstext(sce);
fprintf('\n---- Methods section (copy into your manuscript) ----\n');
fprintf('%s\n', txt);
fprintf('------------------------------------------------------\n\n');
gui.i_showmethodstext(FigureHandle, txt);
end

function i_reportllmword(parentfig, nwritten, failed, folder)
%I_REPORTLLMWORD Tell the user how an LLM Word-report batch went.
%
%   gui.i_reportllmword(parentfig, nwritten, failed, folder) shows one
%   dialog after gui.sc_llm_enrichr2word or gui.sc_llm_dp2word: how many
%   reports were written to FOLDER, and which input files (string array
%   FAILED) produced none. Their reasons are in the Command Window.
%
%   See also gui.sc_llm_enrichr2word, gui.sc_llm_dp2word.

if isempty(failed)
    gui.myHelpdlg(parentfig, sprintf('%s written to %s.', ...
        pkg.i_plural(nwritten, 'Word report'), folder));
    return;
end

msg = sprintf('%s written to %s. No report for: %s. See the Command Window for why.', ...
    pkg.i_plural(nwritten, 'Word report'), folder, strjoin(failed, ', '));
if nwritten == 0
    gui.myWarndlg(parentfig, msg);
else
    gui.myHelpdlg(parentfig, msg);
end
end

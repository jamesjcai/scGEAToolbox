function txt = i_methodsshow(sce, step)
%I_METHODSSHOW Print the Methods paragraph a step stored on SCE.
%
%   pkg.i_methodsshow(sce, step)
%   txt = pkg.i_methodsshow(sce, step) returns the paragraph followed by
%   its references instead of printing, for a GUI to show
%   (GUI.I_SHOWMETHODSTEXT); "" if the step has no text.
%
%   For entry points - GUI callbacks, CLI commands - that run an SCE method
%   which records its paragraph silently (PKG.I_METHODSREPORT with
%   Print=false), because the method is also called inside other analyses
%   where printing would be noise. Prints nothing if the step has no text.
%
%   See also PKG.I_METHODSREPORT, SC_METHODSTEXT.

arguments
    sce
    step (1,1) string
end
txt = "";
if ~isa(sce, 'SingleCellExperiment')
    return
end
stored = sce.struct_methods_text;
group = pkg.i_methodsgroup(string(fieldnames(stored)), step);
if isempty(group), return; end
text = strings(1, numel(group));
refkeys = strings(1, 0);
for j = 1:numel(group)
    text(j) = stored.(group(j)).Paragraph;
    refkeys = [refkeys, stored.(group(j)).References]; %#ok<AGROW>
end
if nargout > 0
    txt = strjoin(text, " ");
    refs = pkg.i_methodsrefs(refkeys);
    if ~isempty(refs)
        txt = txt + string(newline) + newline + "References" + newline + ...
            strjoin(refs, newline);
    end
    return
end
pkg.i_methodsreport([], step, strjoin(text, " "), refkeys);
end

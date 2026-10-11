function i_methodsreport(sce, step, paragraph, refkeys, options)
%I_METHODSREPORT Print a step's Methods paragraph and keep it with the data.
%
%   pkg.i_methodsreport(sce, step, paragraph, refkeys)
%   pkg.i_methodsreport(sce, step, paragraph, refkeys, Print=false)
%
%   Prints PARAGRAPH and the references named by REFKEYS
%   (PKG.I_METHODSREFS) between markers, for the user to paste into a
%   manuscript, and stores both in SCE.struct_methods_text.(STEP), where
%   SC_METHODSTEXT collects every step into one Methods section. Pass []
%   for SCE when there is no object - a matrix-level function - and the
%   text is printed only.
%
%   Call it once per user-level run, from the entry point (an SCE method,
%   a GUI callback, a CLI command), never from a function that other
%   analyses call internally: SC_NORM and RUN.ML_ALONA run inside many
%   steps, and would print every time.
%
%   STEP is one of the keys SC_METHODSTEXT orders the section by. A step
%   run again replaces its earlier paragraph.
%
%   See also SC_METHODSTEXT, PKG.I_METHODSREFS.

arguments
    sce
    step (1,1) string
    paragraph (1,1) string
    refkeys (1,:) string = strings(1, 0)
    options.Print (1,1) logical = true
end

if isa(sce, 'SingleCellExperiment') && isscalar(sce) && isvalid(sce)
    sce.struct_methods_text.(step) = struct('Paragraph', paragraph, ...
        'References', refkeys, 'Time', datetime('now'));
end

if options.Print
    refs = pkg.i_methodsrefs(refkeys);
    fprintf('\n---- Methods text (copy into your manuscript) ----\n');
    fprintf('%s\n', paragraph);
    if ~isempty(refs)
        fprintf('\nReferences\n');
        fprintf('%s\n', refs);
    end
    fprintf('---------------------------------------------------\n\n');
end
end

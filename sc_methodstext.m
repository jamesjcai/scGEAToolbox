function txt = sc_methodstext(sce)
%SC_METHODSTEXT Methods section for the analysis steps run on a dataset.
%
%   sc_methodstext(sce)        prints it, ready to paste into a manuscript
%   txt = sc_methodstext(sce)  returns it instead
%
%   Each analysis step that writes Methods text - quality control, feature
%   selection and embedding, batch correction, clustering, cell type
%   annotation, differential expression, variability and programs, and
%   gene regulatory networks - prints its own paragraph when it runs and
%   keeps a copy in sce.struct_methods_text. This collects those copies
%   in workflow order, each under a run-in heading, followed by one merged
%   reference list. A step run more than once is described by its last
%   run.
%
%   The text says what the toolbox did with the settings that were used;
%   it does not know why. Read it before pasting, and add what only the
%   author knows: the data, the reasons for the thresholds, and any step
%   done outside scGEAToolbox.
%
%   EXAMPLE
%     sce = sce.qcfilter();
%     sce = sce.embedcells('umap');
%     sce = sce.clustercells(10);
%     sce = sce.assigncelltype('mouse');
%     sc_methodstext(sce)
%
%   See also PKG.I_METHODSREPORT, PKG.I_METHODSREFS.

arguments
    sce (1,1) SingleCellExperiment
end

steps = i_steps();
stored = sce.struct_methods_text;
paragraphs = strings(0, 1);
refkeys = strings(1, 0);
names = string(fieldnames(stored));
used = false(size(names));
for k = 1:height(steps)
    % A step may be stored under sub-keys - embedding_features,
    % embedding_umap3d - which share its heading.
    group = pkg.i_methodsgroup(names, steps.Key(k));
    if isempty(group), continue; end
    used = used | ismember(names, group);
    text = strings(1, numel(group));
    for j = 1:numel(group)
        rec = stored.(group(j));
        text(j) = rec.Paragraph;
        refkeys = [refkeys, rec.References]; %#ok<AGROW>
    end
    paragraphs(end+1, 1) = steps.Heading(k) + ". " + strjoin(text, " "); %#ok<AGROW>
end
% A step not in the table above still belongs in the section, after the
% known ones, rather than being dropped silently.
other = names(~used);
for k = 1:numel(other)
    rec = stored.(other(k));
    paragraphs(end+1, 1) = rec.Paragraph; %#ok<AGROW>
    refkeys = [refkeys, rec.References]; %#ok<AGROW>
end

if isempty(paragraphs)
    text = "No Methods text has been recorded for this dataset yet. " + ...
        "Steps record it as they run: quality control, embedding, " + ...
        "clustering, cell type annotation and the analyses.";
else
    % Alphabetical, by first author, as a reference list is read.
    refs = sort(pkg.i_methodsrefs(refkeys));
    % string(newline) first: newline + newline adds two chars as numbers.
    text = strjoin(paragraphs, string(newline) + newline);
    if ~isempty(refs)
        text = text + newline + newline + "References" + newline + ...
            strjoin(refs, newline);
    end
end

if nargout > 0
    txt = text;
else
    fprintf('\n---- Methods section (copy into your manuscript) ----\n');
    fprintf('%s\n', text);
    fprintf('------------------------------------------------------\n\n');
end
end


function steps = i_steps()
% Workflow order, and the run-in heading each step's paragraph gets in the
% combined section.
Key = ["qc"; "embedding"; "harmony"; "clustering"; "annotation"; ...
    "de"; "dv"; "dp"; "grn"; "knk"; "ccc"];
Heading = ["Quality control"; "Feature selection and embedding"; ...
    "Batch correction"; "Clustering"; "Cell type annotation"; ...
    "Differential expression"; "Differential variability"; ...
    "Differential gene programs"; "Gene regulatory networks"; ...
    "Virtual gene knockout"; "Cell-cell communication"];
steps = table(Key, Heading);
end

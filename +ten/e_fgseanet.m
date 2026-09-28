function [OUT] = e_fgseanet(Tf, varargin)
% Merge similar gene sets (Jaccard index > cutoff) in fGSEA report

p = inputParser;
addOptional(p, 'JaccardCutoff', 0.6, @(x) x > 0 & x < 1);
addOptional(p, 'PlotNetwork', true, @islogical);
addOptional(p, 'ShowNotepad', true, @islogical);
parse(p, varargin{:});
jaccardcutoff = p.Results.JaccardCutoff;
plotnetwork = p.Results.PlotNetwork;
shownotepad = p.Results.ShowNotepad;
if size(Tf, 1) < 5
    error('fGSEA output table is too short.')
end

n = size(Tf.leadingEdge, 1);
A = zeros(n);
for i = 1:n - 1
    for j = i + 1:n
        a = strsplit(Tf.leadingEdge{i}, ";");
        b = strsplit(Tf.leadingEdge{j}, ";");
        A(i, j) = length(intersect(a, b)) ./ length(unique(union(a, b)));
        A(j, i) = A(i, j);
    end
end

%%
nodenames = Tf.pathway;
nodenamesfull = Tf.pathway;
for k = 1:n
    nodenamesfull{k} = sprintf('%d_%s', k, Tf.pathway{k});
    nodenames{k} = sprintf('%d', k);
end

%%
B = A .* (A > jaccardcutoff);
G = digraph(B, nodenames);

if plotnetwork
    % gui.myFigure, so a caller can draw it into an open window with
    % gui.myFigure.drawInto.
    hx = gui.myFigure();
    plot(hx.AxHandle, G, 'NodeLabel', nodenames, 'NodeLabelMode', 'auto');
    hx.setTitle('GSEA term network');
    hx.show();
end

%%
[bins, binsizes] = conncomp(G);
[~, idx] = sort(binsizes, 'descend');
OUT = cell(max(bins), 2);


tmpName = [tempname, '.txt'];
fid = fopen(tmpName, 'w');
for k = 1:max(bins)
    fprintf(fid, '\nEnriched Function Group %d\n', k);
    vi = find(bins == idx(k));
    Gx = [];
    for kk = 1:length(vi)
        a = strsplit(Tf.leadingEdge{vi(kk)}, ";");
        Gx = [Gx, a];
    end
    Gx = unique(Gx, 'stable');
    OUT{k, 1} = string(Gx);

    fprintf(fid, '%s ', string(Gx));
    fprintf(fid, '\n');
    fprintf(fid, '\t%s\n', nodenamesfull{bins == idx(k)});
    OUT{k, 2} = deblank(sprintf('%s\n', nodenamesfull{bins == idx(k)}));
end
fclose(fid);

if shownotepad
    [status] = system(['notepad "', tmpName, '" &']);
    if status ~= 0
        if ~(ismcc || isdeployed)
            %#exclude edit
            edit(tmpName);
        end
    end
end

% consider using https://bioconductor.org/packages/release/bioc/vignettes/simplifyEnrichment/inst/doc/simplifyEnrichment.html
end

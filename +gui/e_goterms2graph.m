function e_goterms2graph(filename, varargin)

    if nargin < 1 || isempty(filename)
        [file, path] = uigetfile({'*.xlsx;*.xls', 'Excel Files (*.xlsx, *.xls)'; ...
                                  '*.*', 'All Files (*.*)'}, ...
                                  'Select an Excel File');
        
        if isequal(file, 0)
            error('No file selected.'); % Handle cancel action
        end
        filename = fullfile(path, file);
    end

p = inputParser;
addOptional(p, 'JaccardCutoff', 0.4, @(x) x > 0 & x < 1)
addOptional(p, 'PlotNetwork', true, @islogical);
addOptional(p, 'ShowNotepad', false, @islogical);
parse(p, varargin{:});
jaccardcutoff = p.Results.JaccardCutoff;
plotnetwork = p.Results.PlotNetwork;
shownotepad = p.Results.ShowNotepad;

%cd 'G:\My Drive\Collab\scRNAseq\2503_Xiaofang\test\scgeatool_DEVPAnalysis_Batch_workingfolder'
%filename = 'DEVP_DV_ko2_vs_wt2_Fibroblasts.xlsx';
%if nargin<2, 
%    sheetname = []; % 'Dn_250_Reactome';

%if isempty(sheetname)
    sheets = sheetnames(filename);
%end
    
    % 'Up_250_GO_BP'; end
for kx = 1:length(sheets)
    Tf=readtable(filename,'Sheet', sheets(kx));
    v = Tf.Properties.VariableNames;
    if ~(contains('OverlappingGenes', v) && contains('TermName', v))
        continue;
    end
    
    n = size(Tf.OverlappingGenes, 1);
    A = zeros(n);
    for i = 1:n - 1
        for j = i + 1:n
            a = strsplit(Tf.OverlappingGenes{i}, ";");
            b = strsplit(Tf.OverlappingGenes{j}, ";");
            A(i, j) = length(intersect(a, b)) ./ length(unique(union(a, b)));
            A(j, i) = A(i, j);
        end
    end
    
    if size(Tf, 1) < 5
        error('Table is too short.')
    end
    
    
    %%
    nodenames = Tf.TermName;
    nodenamesfull = Tf.TermName;
    nodenameS = Tf.TermName;
    for k = 1:n
        % nodenamesfull{k}=sprintf('%d_%s',k,Tf.TermName{k});
        % nodenamesfull{k}=sprintf('%s',Tf.TermName{k});
        nodenamesfull{k} = sprintf('%d_%s', k, Tf.TermName{k});
        nodenameS{k} = sprintf('%s', Tf.TermName{k});
        %a=sprintf('%d\\_%s',k,Tf.TermName{k});
        %a=extractBefore(a,min(20,length(a)));
        nodenames{k} = sprintf('%d', k);
    end
    
    %%
    %B=A.*(abs(A)>quantile(abs(A(:)),0.95));
    B = A .* (A > jaccardcutoff);
    % G=digraph(A,Tf.TermName);
    G = graph(B, nodenames);
    % LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
    % LWidths(LWidths==0)=1e-5;
    
    if plotnetwork
        % figure;
        % p = plot(G, 'NodeLabel', nodenames, ...
        %     'NodeLabelMode', 'auto');
    
        % sc_grnview(B, nodenameS, [], []);
        G = pkg.i_makegraph(B, nodenameS);
        gui.i_singlegraph(G, strrep(sheets(kx),'_','\_'));
    end
end

% p=plot(G,'NodeLabel',nodenames,'NodeFontAngle','normal',...
%     'NodeFontSize',12);
% if ~isempty(LWidths)
%     p.LineWidth=LWidths;
% end
% p.MarkerSize = 7;
% p.Marker = 's';
% p.NodeColor = 'r';

%%

if shownotepad
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
        a = strsplit(Tf.OverlappingGenes{vi(kk)}, ";");
        Gx = [Gx, a];
    end
    Gx = unique(Gx, 'stable');
    OUT{k, 1} = string(Gx);

    fprintf(fid, '%s ', string(Gx));
    fprintf(fid, '\n');
    fprintf(fid, '\t%s\n', nodenamesfull{bins == idx(k)});
    OUT{k, 2} = deblank(sprintf('%s\n', nodenamesfull{bins == idx(k)}));
end
%fprintf(fid,'---------------\n');
fclose(fid);


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

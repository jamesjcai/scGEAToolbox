function [requirerefresh] = callback_GenerateSyntheticCells(src, ptsSelected)

requirerefresh = false;
if nargin < 1, ptsSelected = []; end
if isempty(ptsSelected), return; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

nBrushed = nnz(ptsSelected);
nTotal = sce.NumCells;  % gene-by-cell matrix assumption

msg = sprintf(['Generate %d synthetic cells based on %d brushed cells using generative Fourier transformer [PMID: 39843603]? ' ...
               'Dataset size will become %d cells.'], ...
               nBrushed, nBrushed, nTotal + nBrushed);

answer = gui.myQuestdlg(FigureHandle, msg);
if ~strcmp(answer,'Yes'), return; end



tag = 'is_synthetic';
issynthetic = sce.getCellAttribute(tag);
if isempty(issynthetic)
    issynthetic = false(sce.NumCells, 1);
    sce.setCellAttribute(tag, issynthetic);
else
    if any(issynthetic(ptsSelected))
        % warning('Synthetic cells are selected.');
    end
end

oldcn = sce.NumCells;
oldgn = sce.NumGenes;

if sce.NumCells*sce.NumGenes < 4e8
    sceori = copy(sce);
else
    answer = gui.myQuestdlg(FigureHandle, ...
        'You are about to change the SCE data. This cannot be undone.');
    if ~strcmp(answer, 'Yes'), return; end        
    sceori = [];
end


%fw = gui.myWaitbar(FigureHandle);

sub_sce = copy(sce).selectcells(ptsSelected);
x = sce.X(:, ptsSelected);
x = sc_ifft(x);
sub_sce.X = x;
s = sub_sce.s;
r = trnd(100, size(s)) * 0.75;
sub_sce.s = s + r;  % * (max(range(s)));

sub_sce.setCellAttribute(tag, true(sub_sce.NumCells, 1));

sce = sc_mergesces({sce, sub_sce}, "intersect", true, true);

f = fieldnames(sce.struct_cell_embeddings);

for i = 1:numel(f)
    if ~isempty(sce.struct_cell_embeddings.(f{i})) 
        s = sce.struct_cell_embeddings.(f{i});
        r = trnd(100, size(s)) * 0.75;
        sce.struct_cell_embeddings.(f{i}) = s + r;
    end
end

%gui.myWaitbar(FigureHandle, fw);

newcn = sce.NumCells;
newgn = sce.NumGenes;

    if oldcn-newcn==0 && oldgn-newgn==0
        gui.myHelpdlg(FigureHandle, "No cells are synthesized.");
        % requirerefresh = false;
        return;
    end
    if ~isempty(sceori)
        answer = gui.myQuestdlg(FigureHandle, ...
            sprintf('%d cells will be added.\n[%d genes x %d cells] => [%d genes x %d cells]', ...
                newcn - oldcn, oldgn, oldcn, newgn, newcn),'', ...
                {'Accept Changes', 'Cancel Changes'}, 'Accept Changes');
        if ~strcmp(answer, 'Accept Changes')
            % requirerefresh = false;
            % sce = sceori;
            sce = copy(sceori);
        else
            requirerefresh = true;
        end
    else
        requirerefresh = true;
    end


    gui.myGuidata(FigureHandle, sce, src);

end

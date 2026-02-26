function [requirerefresh] = callback_GenerateSyntheticCells(src, ptsSelected)

requirerefresh = false;

if nargin < 1, ptsSelected = []; end
if isempty(ptsSelected), return; end

% assignin("base", 'ptsSelected', ptsSelected)

[FigureHandle, sce] = gui.gui_getfigsce(src);

oldcn = sce.NumCells;
oldgn = sce.NumGenes;

if sce.NumCells*sce.NumGenes < 4e8
    sceori = copy(sce);
    % disp('Ready for reversible.');
else
    answer = gui.myQuestdlg(FigureHandle, ...
        'You are about to change the SCE data. This cannot be undone.');
    if ~strcmp(answer, 'Yes'), return; end        
    sceori = [];
end


fw = gui.myWaitbar(FigureHandle);

sub_sce = copy(sce).selectcells(ptsSelected);
x = sce.X(:, ptsSelected);
x = sc_ifft(x);
sub_sce.X = x;
s = sub_sce.s;
sub_sce.s = s + trnd(100, size(s)) * 0.75; % * (max(range(s)));
sce = sc_mergesces({sce, sub_sce}, "intersect", true, true);

gui.myWaitbar(FigureHandle, fw);

newcn = sce.NumCells;
newgn = sce.NumGenes;

    if oldcn-newcn==0 && oldgn-newgn==0
        gui.myHelpdlg(FigureHandle, "No cells are synthesized.");
        % requirerefresh = false;
        return;
    end
    if ~isempty(sceori)
        answer = gui.myQuestdlg(FigureHandle, sprintf('%d cells will be added.\n[%d genes x %d cells] => [%d genes x %d cells]', ...
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

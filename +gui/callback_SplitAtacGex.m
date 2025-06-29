function [sceg, scep] = callback_SplitAtacGex(src, ~)



    [FigureHandle, sce] = gui.gui_getfigsce(src);

if ~any(contains(sce.g, ':'))
    gui.myWarndlg(FigureHandle, 'Not a multiome ATAC+GEX matrix.');
    return;
end

ispeak = contains(sce.g, ':');
isgene = ~ispeak;

sceg = sce;
sceg.g = sceg.g(isgene);
sceg.X = sceg.X(isgene, :);

scep = sce;
scep.g = scep.g(~isgene);
scep.X = scep.X(~isgene, :);


labels = {'Save SCEG to variable named:', ...
    'Save SCEP to variable named:'};
vars = {'sceg', 'scep'};
values = {sceg, scep};
[~, ~] = export2wsdlg(labels, vars, values, ...
    'Save Data to Workspace');

% --------------------
pause(2);
fx = scgeatool(sceg);
fx.Position(3:4) = 0.8 * fx.Position(3:4);
movegui(fx, 'center');
fx.Position(1) = fx.Position(1) - 250;
fy = fx.CurrentAxes;
fy.Subtitle.String = '[genes x cells]';
answer = gui.myQuestdlg(FigureHandle, ...
    'scRNAseq data extracted. Continue?', '');
if isempty(answer), return; end
if ~strcmp(answer, 'Yes'), return; end
fy = scgeatool(scep);
fy.Position(3:4) = 0.8 * fy.Position(3:4);
movegui(fy, 'center');
fy.Position(1) = fy.Position(1) + 250;
fy = fy.CurrentAxes;
fy.Subtitle.String = '[peaks x cells]';

end
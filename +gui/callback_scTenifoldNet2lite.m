function callback_scTenifoldNet2lite(src, ~)

import ten.*
import pkg.*

FigureHandle = src.Parent.Parent;
if ~gui.gui_showrefinfo('scTenifoldNet [PMID:33336197]'), return; end

sce = guidata(FigureHandle);

[i1, i2] = gui.i_select2smplgrps(sce, false);
if isscalar(i1) || isscalar(i2), return; end

fw = gui.gui_waitbar;
disp('Constructing networks (1/2) ...')
X = sc_norm(sce.X);
X = log1p(X);
X0 = X(:, i1);
X1 = X(:, i2);
A0 = sc_pcnetpar(X0);
disp('Constructing networks (2/2) ...')
A1 = sc_pcnetpar(X1);
A0sym = 0.5 * (A0 + A0');
A1sym = 0.5 * (A1 + A1');

disp('Manifold alignment...')
[aln0, aln1] = i_ma(A0sym, A1sym);
disp('Differential regulation (DR) detection...')
glist = sce.g;
T = i_dr(aln0, aln1, glist);
gui.gui_waitbar(fw);

tstr = matlab.lang.makeValidName(string(datetime));
b = 'sctenifoldnet_outs';
a = sprintf('output_%s', tstr);
if ~exist(fullfile(tempdir, b), 'dir')
    mkdir(fullfile(tempdir, b));
    pause(1);
end
f1 = fullfile(tempdir, b, a);
save(f1, 'T');
a = sprintf('output_%s.xlsx', tstr);
f1 = fullfile(tempdir, b, a);
writetable(T, f1, 'FileType', 'spreadsheet');
fprintf('The result has been saved in %s\n', f1);
bx = sprintf('The result has been saved in %s. Open the folder to locate it?', f1);
answer = questdlg(bx);
switch answer
    case 'Yes'
        winopen(fullfile(tempdir, b));
end


%{
figure;
ten.e_mkqqplot(T);
% answer223=questdlg('Run GSEA analysis?');
answer223=gui.i_questdlgtimer(15,'Run GSEA analysis?');
if ~isempty(answer223) && strcmp(answer223,'Yes')
    gseaok=true;
    try
        Tr=ten.e_fgsearun(T);
        save(sprintf('T_GSEAres_%s',tstr),'Tr');
    catch ME
        warning(ME.message);
        gseaok=false;
    end
    if gseaok
        answer323=gui.i_questdlgtimer(15,'Group GSEA hits?');
        if ~isempty(answer323) && strcmp(answer323,'Yes')
            ten.e_fgseanet(Tr);
        end
    end
end
gui.i_exporttable(T,true,'T_DRgenes');
if gseaok
    gui.i_exporttable(Tr,true,'T_GSEAres');
end
%}

end

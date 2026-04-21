function callback_scTenifoldNet1lite(src, events)

import pkg.*

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~gui.gui_showrefinfo('scTenifoldNet [PMID:33336197]', FigureHandle), return; end

extprogname = 'ml_scTenifoldNet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

if numel(unique(sce.c_cell_type_tx)) > 1
    answer = gui.myQuestdlg(FigureHandle, ...
        'Construct gene regulatory network (GRN) for all cells or selected cells?', ...
        '', {'All Cells', 'Select Cells...', 'Cancel'}, 'All Cells');
    switch answer
        case 'Cancel'
            return;
        case 'All Cells'
        case 'Select Cells...'
            gui.callback_SelectCellsByClass(src, events);
            return;
        otherwise
            return;
    end
end

X = sc_norm(sce.X);
X = log1p(X);
useGPU = pkg.i_usegpu(X);
if useGPU
    disp('GPU detected — using CUDA GPU acceleration.');
    useparallel = false;
else
    answer = gui.myQuestdlg(FigureHandle, 'Use parallel computing?', 'Parallel Computing', ...
        {'Use parallel', 'Not use parallel'}, 'Use parallel');
    switch answer
        case 'Use parallel'
            useparallel = true;
        case 'Not use parallel'
            useparallel = false;
        otherwise
            return;
    end
end

fw = gui.myWaitbar(FigureHandle);
disp('Constructing gene regulatory network...')
try
    A = net.pcrnet(X, 3, true, true, useparallel, ~useparallel, useGPU);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

tstr = matlab.lang.makeValidName(string(datetime));
b = 'sctenifoldnet_outs';
a = sprintf('output_%s', tstr);
if ~exist(fullfile(wkdir, b), 'dir')
    mkdir(fullfile(wkdir, b));
    drawnow;
end
f1 = fullfile(wkdir, b, a);
g = sce.g;
save(f1, 'A', 'g', '-v7.3');
fprintf('The network has been saved in %s.mat\n', f1);
bx = sprintf('The network has been saved in %s.mat. Open the folder to locate it?', f1);
if strcmp('Yes', gui.myQuestdlg(FigureHandle, bx))
    winopen(fullfile(wkdir, b));
end

if ~(ismcc || isdeployed)
    labels = {'Save network to variable named:', 'Save sce.g to variable named:'};
    vars = {'A', 'g'};
    values = {A, sce.g};
    export2wsdlg(labels, vars, values);
end

end

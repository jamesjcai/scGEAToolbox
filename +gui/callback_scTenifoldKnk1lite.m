function callback_scTenifoldKnk1lite(src, ~)

import pkg.*

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('scTenifoldKnk [PMID:35510185]', FigureHandle), return; end

extprogname = 'scTenifoldKnk';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

answer = gui.myQuestdlg(FigureHandle, ...
    'Construct network de novo or use existing network in Workspace?', ...
    'Input Network', {'Construct de novo', 'Use existing'}, 'Construct de novo');
switch answer
    case 'Use existing'
        A0 = in_loadnetwork(FigureHandle, sce);
        if isempty(A0), return; end
    case 'Construct de novo'
        A0 = [];
        gui.myHelpdlg(FigureHandle, "Network will be constructed. Now, " + ...
            "select a KO gene (i.e., gene to be knocked out).");
    otherwise
        return;
end

gsorted = natsort(sce.g);
if isempty(gsorted), return; end

if gui.i_isuifig(FigureHandle)
    [indx2, tf] = gui.myListdlg(FigureHandle, gsorted, 'Select a KO gene');
else
    [indx2, tf] = listdlg('PromptString', {'Select a KO gene'}, ...
        'SelectionMode', 'single', 'ListString', gsorted, 'ListSize', [220, 300]);
end
if tf ~= 1, return; end
[~, idx] = ismember(gsorted(indx2), sce.g);

if isempty(A0)
    answer2 = gui.myQuestdlg(FigureHandle, ...
        sprintf('Ready to construct network and then knock out %s (gene #%d). Continue?', ...
        sce.g(idx), idx));
    if ~strcmpi(answer2, 'Yes'), return; end

    X = sc_norm(sce.X);
    X = log1p(X);
    useGPU = pkg.i_usegpu(X);
    if useGPU
        disp('GPU detected — using CUDA GPU acceleration.');
        useparallel = false;
    else
        answer3 = gui.myQuestdlg(FigureHandle, 'Use parallel computing?', 'Parallel Computing', ...
            {'Use parallel', 'Not use parallel'}, 'Use parallel');
        switch answer3
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
        A0 = net.pcrnet(X, 3, true, true, useparallel, ~useparallel, useGPU);
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    gui.myWaitbar(FigureHandle, fw);

    if ~(ismcc || isdeployed)
        labels = {'Save network to variable named:'};
        vars = {'A0'};
        values = {A0};
        export2wsdlg(labels, vars, values);
    end
else
    answer2 = gui.myQuestdlg(FigureHandle, ...
        sprintf('Ready to knock out %s (gene #%d). Continue?', sce.g(idx), idx));
    if ~strcmpi(answer2, 'Yes'), return; end
end

if nnz(A0(idx, :) ~= 0) == 0
    gui.myWarndlg(FigureHandle, sprintf('KO gene (%s) has no links with other genes.', sce.g(idx)));
    return;
elseif nnz(A0(idx, :) ~= 0) < 50
    s = sprintf('KO gene (%s) has too few links (n=%d). Continue?', ...
        sce.g(idx), nnz(A0(idx, :) ~= 0));
    if ~strcmpi(gui.myQuestdlg(FigureHandle, s, '', [], [], 'error'), 'Yes'), return; end
end

fw = gui.myWaitbar(FigureHandle);
disp('>> [T] = ten.i_knk(A0, idx, sce.g, true);')
try
    T = ten.i_knk(A0, idx, sce.g, true);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

[answer_exp, filename] = gui.i_exporttable(T, true, ...
    sprintf('Ttenifldknk_%s', sce.g(idx)), ...
    sprintf('TenifldKnkTable_%s', sce.g(idx)), [], [], FigureHandle);
if isempty(filename)
    fprintf('\nResults have been saved in %s.\n\n', answer_exp);
else
    fprintf('\nResults have been saved in %s: %s.\n\n', answer_exp, filename);
end
disp('Downstream Analysis Options:');
disp('===============================');
disp('run.web_Enrichr(T.genelist(1:200));');
disp('Tf=ten.e_fgsearun(T);');
disp('Tn=ten.e_fgseanet(Tf);');
disp('===============================');


    function A0 = in_loadnetwork(hfig, sce_)
        A0 = [];
        a = evalin('base', 'whos');
        valididx = false(length(a), 1);
        for k = 1:length(a)
            if max(a(k).size) == sce_.NumGenes && min(a(k).size) == sce_.NumGenes
                valididx(k) = true;
            end
        end
        if isempty(a) || ~any(valididx)
            if ~strcmpi(gui.myQuestdlg(hfig, ...
                    'Workspace contains no network variable. Read from .mat file?', ''), 'Yes')
                return;
            end
            A0 = in_readfile(hfig, sce_.NumGenes);
        else
            a = a(valididx);
            b = struct2cell(a);
            if gui.i_isuifig(hfig)
                [indx, tf] = gui.myListdlg(hfig, b(1, :), 'Select network variable:');
            else
                [indx, tf] = listdlg('PromptString', {'Select network variable:'}, ...
                    'liststring', b(1, :), 'SelectionMode', 'single', 'ListSize', [220, 300]);
            end
            if tf ~= 1, return; end
            A0 = evalin('base', a(indx).name);
            if size(A0, 1) ~= size(A0, 2) || size(A0, 2) ~= sce_.NumGenes
                gui.myErrordlg(hfig, 'Not a valid network.');
                A0 = [];
            end
        end
    end

    function A0 = in_readfile(hfig, n)
        A0 = [];
        if gui.i_isuifig(hfig)
            [fname, pathname] = uigetfile(hfig, ...
                {'*.mat', 'Saved GRN Files (*.mat)'; '*.*', 'All Files (*.*)'}, ...
                'Pick a GRN Data File');
        else
            [fname, pathname] = uigetfile( ...
                {'*.mat', 'Saved GRN Files (*.mat)'; '*.*', 'All Files (*.*)'}, ...
                'Pick a GRN Data File');
        end
        if isequal(fname, 0), return; end
        data = load(fullfile(pathname, fname), 'A0');
        try
            A0 = data.A0;
        catch ME
            disp(ME.message);
            return;
        end
        if size(A0, 1) ~= n || size(A0, 2) ~= n
            A0 = [];
        end
    end

end

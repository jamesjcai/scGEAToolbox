function callback_scTenifoldKnk1(src, ~)

[FigureHandle] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('scTenifoldKnk [PMID:35510185]', ...
        FigureHandle), return; end

try
    ten.check_tensor_toolbox;
catch
    gui.i_installtensortoolbox(src);
    try
        ten.check_tensor_toolbox;
    catch ME
        gui.myErrordlg(FigureHandle, ME.message);
        return;
    end
end

extprogname = 'scTenifoldKnk';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
olddir = pwd;
cleanupObj = onCleanup(@() cd(olddir));
if isfolder(wkdir), cd(wkdir); end
import ten.*

[FigureHandle, sce] = gui.gui_getfigsce(src);
netname = '';   % where an existing network came from, for the confirm prompt

answer = gui.myQuestdlg(FigureHandle, 'Construct network de novo or use existing network in Workspace?', ...
    'Input Network', {'Construct de novo', 'Use existing'}, 'Construct de novo');
switch answer
    case 'Use existing'
        a = evalin('base', 'whos');
        b = struct2cell(a);
        valididx = false(length(a), 1);
        for k = 1:length(a)
            if max(a(k).size) == sce.NumGenes && min(a(k).size) == sce.NumGenes
                valididx(k) = true;
            end
        end
        if isempty(b) || ~any(valididx)
            [anw] = gui.myQuestdlg(FigureHandle, 'Workspace contains no network variable. Read from .mat file?','');
            if ~strcmp(anw, 'Yes'), return; end
            [A0] = in_readA0fromfile(sce.NumGenes);
            netname = 'file';
            if isempty(A0) || size(A0, 1) ~= sce.NumGenes || size(A0, 2) ~= sce.NumGenes
                gui.myErrordlg(FigureHandle, 'Not a valid network.');
                return;
            end
        else
            a = a(valididx);
            b = b(:, valididx);

       if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, b(1, :), ...
                'Select network variable:', [], false);
        else
            [indx, tf] = listdlg('PromptString', {'Select network variable:'}, ...
                'liststring', b(1, :), 'SelectionMode', 'single', 'ListSize', [220, 300]);
        end

            if tf == 1
                A0 = evalin('base', a(indx).name);
                netname = a(indx).name;
            else
                return;
            end
            [m, n] = size(A0);
            if m ~= n || n ~= length(sce.g)
                gui.myErrordlg(FigureHandle, 'Not a valid network.');
                return;
            end
        end
    case 'Construct de novo'
        try
            ten.check_tensor_toolbox;
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        A0 = [];
        gui.myHelpdlg(FigureHandle, "Network will be constructed. Now, " + ...
            "select a KO gene (i.e., gene to be knocked out).");
    otherwise
        return;
end
gsorted = natsort(sce.g);
if isempty(gsorted), return; end

   if gui.i_isuifig(FigureHandle)
        [indx2, tf] = gui.myListdlg(FigureHandle, gsorted, 'Select a KO gene', [], false);
    else
        [indx2, tf] = listdlg('PromptString', {'Select a KO gene'}, ...
            'SelectionMode', 'single', 'ListString', gsorted, 'ListSize', [220, 300]);
    end

if tf == 1
    [~, idx] = ismember(gsorted(indx2), sce.g);
else
    return;
end


if isempty(A0)
    answer = gui.myQuestdlg(FigureHandle, sprintf('Ready to construct network and then knock out %s (gene #%d). Continue?', ...
        sce.g(idx), idx));
else
    answer = gui.myQuestdlg(FigureHandle, sprintf('Ready to knock out %s (gene #%d) from network (%s). Continue?', ...
        sce.g(idx), idx, netname));
end

if ~strcmpi(answer, 'Yes'), return; end

if isempty(A0)
    [nsubsmpl, csubsmpl, savegrn] = gui.i_tenifoldnetpara(FigureHandle);
    if isempty(nsubsmpl) || isempty(csubsmpl) || isempty(savegrn), return; end
    try
        fw = gui.myWaitbar(FigureHandle);
        parfor k=1:32
        end
        if pkg.i_isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end

        [T, A0] = ten.sctenifoldknk(sce.X, sce.g, idx, ...
            'sorttable', true, 'nsubsmpl', nsubsmpl, 'csubsmpl', csubsmpl, ...
            'savegrn', savegrn);
        gui.myWaitbar(FigureHandle, fw);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    isreconstructed = true;
else
    if nnz(A0(:, idx) ~= 0) == 0
        s = sprintf('KO gene (%s) has no link or too few links (n<50) with other genes.', ...
            sce.g(idx));
        gui.myWarndlg(FigureHandle, s);
        return;
    elseif nnz(A0(:, idx) ~= 0) < 50
        s = sprintf('KO gene (%s) has too few links (n=%d) with other genes. Continue?', ...
            sce.g(idx), nnz(A0(:, idx) ~= 0));
        answer11 = gui.myQuestdlg(FigureHandle, s,'',[],[],'error');
        switch answer11
            case 'Yes'
                doit = true;
            case 'No'
                return;
            case 'Cancel'
                return;
            otherwise
                return;
        end
    else
        doit = true;
    end

    if doit
        try
            fw = gui.myWaitbar(FigureHandle);
            disp('>> [T] = ten.i_knk(A0, targetgene, genelist, true);')
            [T] = ten.i_knk(A0, idx, sce.g, true);
            gui.myWaitbar(FigureHandle, fw);
        catch ME
            gui.myWaitbar(FigureHandle, fw);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
    end
    isreconstructed = false;
end

% Results first, then the network: the table is what was asked for.
[answer, filename] = gui.i_exporttable(T, true, ...
    sprintf('Ttenifldknk_%s', sce.g(idx)), ...
    sprintf('TenifldKnkTable_%s', sce.g(idx)), [], [], FigureHandle);
if isempty(filename)
    fprintf('\nResults have been saved in %s.\n\n', answer);
else
    fprintf('\nResults have been saved in %s: %s.\n\n', answer, filename);
end
disp('Downstream Analysis Options:');
disp('===============================');
disp('run.web_Enrichr(T.genelist(1:200));');
disp('Tf=ten.e_fgsearun(T);');
disp('Tn=ten.e_fgseanet(Tf);');
disp('===============================');

if isreconstructed
    in_offersavenetwork(A0);
end

    function [A0] = in_readA0fromfile(n)
        A0 = [];
        % uigetfile takes no parent figure (a figure passed first is read
        % as the filter spec); raise the app afterwards.
        [fname, pathname] = uigetfile( ...
            {'*.mat', 'Saved GRN Files (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Pick a GRN Data File');
        if pkg.i_isvalid(FigureHandle), figure(FigureHandle); end
        if isequal(fname, 0), return; end
            filen = fullfile(pathname, fname);
            data = load(filen, 'A0');

            try
                A0 = data.A0;
            catch ME
                disp(ME.message);
            end
            if ~isempty(A0)
                if ~(size(A0,1)==n && size(A0,2) == n)
                    A0 = [];
                end
            end
    end

    function in_offersavenetwork(net)
        % Last step, after the results: every progress bar is closed by now,
        % and the dialog blocks, so nothing else opens on top of it.
        if ismcc || isdeployed, return; end
        labels = {'Save constructed network to variable named:'};
        if gui.i_isuifig(FigureHandle)
            gui.myExport2wsdlg(labels, {'A0'}, {net}, ...
                'Save Network to Workspace', [], FigureHandle);
        else
            waitfor(export2wsdlg(labels, {'A0'}, {net}, ...
                'Save Network to Workspace'));
        end
    end

end

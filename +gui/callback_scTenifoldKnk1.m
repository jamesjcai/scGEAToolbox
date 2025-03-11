function callback_scTenifoldKnk1(src, ~)

if ~gui.gui_showrefinfo('scTenifoldKnk [PMID:35510185]'), return; end

try
    ten.check_tensor_toolbox;
catch ME
    errordlg(ME.message);
    return;
end

extprogname = 'scTenifoldKnk';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end
if isfolder(wkdir)
    cd(wkdir);
end
import ten.*

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, 'Construct network de novo or use existing network in Workspace?', ...
    'Input Network', {'Construct de novo', 'Use existing'}, 'Construct de novo');
switch answer
    case 'Use existing'
        a = evalin('base', 'whos');
        b = struct2cell(a);
        % if isempty(b)
        %     gui.myHelpdlg(FigureHandle, 'No variable in the WorkSpace.', '');
        % 
        %     return;
        % end
        valididx = false(length(a), 1);
        for k = 1:length(a)
            if max(a(k).size) == sce.NumGenes && min(a(k).size) == sce.NumGenes
                valididx(k) = true;
            end
        end
        if isempty(b) || ~any(valididx)
            [anw] = gui.myQuestdlg(FigureHandle, 'Workspace contains no network varible. Read from .mat file?','');
            if ~strcmp(anw, 'Yes'), return; end
            [A0] = in_readA0fromfile(sce.NumGenes);
            if isempty(A0) || size(A0, 1) ~= sce.NumGenes || size(A0, 2) ~= sce.NumGenes
                errordlg('Not a valid network.');
                return;
            end
        else
            %valididx=ismember(b(4,:),'double');
            a = a(valididx);
            b = b(:, valididx);
            [indx, tf] = listdlg('PromptString', {'Select network variable:'}, ...
                'liststring', b(1, :), 'SelectionMode', 'single', 'ListSize', [220, 300]);
            if tf == 1
                A0 = evalin('base', a(indx).name);
            else
                return;
            end
            [m, n] = size(A0);
            if m ~= n || n ~= length(sce.g)
                errordlg('Not a valid network.');
                return;
            end
        end
    case 'Construct de novo'
        try
            ten.check_tensor_toolbox;
        catch ME
            errordlg(ME.message);
            return;
        end
        A0 = [];
        gui.myHelpdlg(FigureHandle, "Network will be constructed. Now, select a KO gene (i.e., gene to be knocked out).", '');
    otherwise
        return;
end
gsorted = natsort(sce.g);
if isempty(gsorted), return; end
[indx2, tf] = listdlg('PromptString', {'Select a KO gene'}, ...
    'SelectionMode', 'single', 'ListString', gsorted, 'ListSize', [220, 300]);
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
        sce.g(idx), idx, a(indx).name));
end

if ~strcmpi(answer, 'Yes'), return; end


if isempty(A0)
    try
        fw = gui.gui_waitbar;
        parfor k=1:32
        end
        figure(FigureHandle);
 
        [T, A0] = ten.sctenifoldknk(sce.X, sce.g, idx, ...
            'sorttable', true, 'nsubsmpl', 10, ...
            'savegrn', isfolder(wkdir));
        gui.gui_waitbar(fw);
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    isreconstructed = true;
else
    doit = false;
    if nnz(A0(idx, :) ~= 0) == 0
        s = sprintf('KO gene (%s) has no link or too few links (n<50) with other genes.', ...
            sce.g(idx));
        gui.myWarndlg(FigureHandle, s);
        return;
    elseif nnz(A0(idx, :) ~= 0) < 50
        s = sprintf('KO gene (%s) has too few links (n=%d) with other genes. Continue?', ...
            sce.g(idx), nnz(A0(idx, :) ~= 0));
        answer11 = gui.myQuestdlg(FigureHandle, s);
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
            fw = gui.gui_waitbar;
            disp('>> [T] = ten.i_knk(A0, targetgene, genelist, true);')
            [T] = ten.i_knk(A0, idx, sce.g, true);
            gui.gui_waitbar(fw);
        catch ME
            gui.gui_waitbar(fw);
            errordlg(ME.message);
            return;
        end
    end
    %A1=A0;
    %A1(idx,:)=0;
    %[aln0,aln1]=i_ma(A0,A1);
    %T=i_dr(aln0,aln1,sce.g,true);
    isreconstructed = false;
end

if ~(ismcc || isdeployed)
    if isreconstructed
        labels = {'Save network to variable named:'};
        vars = {'A0'};
        values = {A0};
        waitfor(export2wsdlg(labels, vars, values));
    end
end


[answer, filename] = gui.i_exporttable(T, true, ...
    sprintf('Ttenifldknk_%s', sce.g(idx)), ...
    sprintf('TenifldKnkTable_%s', sce.g(idx)));
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

    function [A0] = in_readA0fromfile(n)
        A0 = [];
        [fname, pathname] = uigetfile( ...
            {'*.mat', 'Saved GRN Files (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Pick a GRN Data File');
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

end


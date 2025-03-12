function callback_scTenifoldNet2(src, ~)

import ten.*
import pkg.*
try
    ten.check_tensor_toolbox;
catch ME
    errordlg(ME.message);
    return;
end

%     if exist('sctenifoldnet','file')~=2
%         errordlg('scTenifoldNet is not installed.');
%         disp('To install scTenifoldNet, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
%         disp('addpath(''./scTenifoldNet-master/MATLAB'');');
%         return;
%     end

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('scTenifoldNet [PMID:33336197]', FigureHandle), return; end
    


extprogname = 'scTenifoldNet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
if isfolder(wkdir)
    cd(wkdir);
end


answer = gui.myQuestdlg(FigureHandle, 'Construct networks de novo or use existing networks in Workspace?', ...
    'Input Networks', {'Construct de novo', 'Use existing'}, 'Construct de novo');
switch answer
    case 'Use existing'
        a = evalin('base', 'whos');
        b = struct2cell(a);
        %valididx=ismember(b(4,:),'double');
        %a=a(valididx);
        if isempty(b)
            gui.myHelpdlg(FigureHandle, 'No variable in the WorkSpace.');
            return;
        end
        [indx, tf] = listdlg('PromptString', {'Select two networks:'}, ...
            'liststring', b(1, :), ...
            'SelectionMode', 'multiple', ...
            'ListSize', [220, 300]);
        if tf ~= 1, return; end
        if length(indx) ~= 2, gui.myWarndlg(FigureHandle, 'Need two networks.');
            return;
        end
        A0 = evalin('base', a(indx(1)).name);
        A1 = evalin('base', a(indx(2)).name);
        if ~isequal(size(A0), size(A1))
            errordlg('Two networks should be in the same size.');
            return;
        end
        [indx, tf] = listdlg('PromptString', {'Select the gene list:'}, ...
            'liststring', b(1, :), ...
            'SelectionMode', 'single', ...
            'ListSize', [220, 300]);
        if tf ~= 1, return; end
        glist = evalin('base', a(indx).name);
        if length(glist) ~= size(A0, 1)
            errordlg('Not a valid gene list.');
            return;
        end
        fw = gui.myWaitbar(FigureHandle);
        disp('Reading networks...')
        A0sym = 0.5 * (A0 + A0');
        A1sym = 0.5 * (A1 + A1');

        disp('Manifold alignment...')
        [aln0, aln1] = i_ma(A0sym, A1sym);
        disp('Differential regulation (DR) detection...')
        T = i_dr(aln0, aln1, glist);
        gui.myWaitbar(FigureHandle, fw);

    case 'Construct de novo'
        [i1, i2] = gui.i_select2smplgrps(sce, false, FigureHandle);
        if isscalar(i1) || isscalar(i2), return; end

        [nsubsmpl, csubsmpl, savegrn] = gui.i_tenifoldnetpara;
        if isempty(nsubsmpl) || isempty(csubsmpl) || isempty(savegrn), return; end
        if csubsmpl >= min([size(sce.X(:, i1), 2), size(sce.X(:, i2), 2)])
            errordlg('csubsmpl should be smaller than sce.NumCells.');
            return;
        end

        answer123 = gui.myQuestdlg(FigureHandle, 'This analysis may take several hours. Continue?');
        if ~strcmpi(answer123, 'Yes'), return; end

        fw = gui.myWaitbar(FigureHandle);
        try
            fprintf('\n');
            % disp('[T]=ten.sctenifoldnet(X1,X2,g,''nsubsmpl'',10,''csubsmpl'',500,''savegrn'',true);')
            disp('[T]=ten.sctenifoldnet(sce.X(:,idx1),sce.X(:,idx2),sce.g,''nsubsmpl'',10,''csubsmpl'',500,''savegrn'',true);')
            [T] = ten.sctenifoldnet(sce.X(:, i1), sce.X(:, i2), sce.g, ...
                'nsubsmpl', nsubsmpl, 'csubsmpl', csubsmpl, 'savegrn', savegrn);
        catch ME
            gui.myWaitbar(FigureHandle, fw);
            errordlg(ME.message);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
    otherwise
        return;
end


tstr = matlab.lang.makeValidName(string(datetime));
save(sprintf('T_DRgenes_%s', tstr), 'T');
fprintf('The result has been saved in T_DRgenes_%s.mat\n', tstr);

%figure;
ten.e_mkqqplot(T);
answer223=gui.myQuestdlg(FigureHandle, 'Run GSEA analysis?');
if ~isempty(answer223) && strcmp(answer223, 'Yes')
    gseaok = true;
    try
        Tr = ten.e_fgsearun(T);
        save(sprintf('T_GSEAres_%s', tstr), 'Tr');
    catch ME
        warning(ME.message);
        gseaok = false;
    end
    if gseaok
        % answer323 = gui.i_questdlgtimer(15, 'Group GSEA hits?');
        answer323 = gui.myQuestdlg(FigureHandle, 'Group GSEA hits?');
        if ~isempty(answer323) && strcmp(answer323, 'Yes') && ~isempty(Tr)
            ten.e_fgseanet(Tr);
        end
    end
end
gui.i_exporttable(T, true, 'Ttenifldnet', 'TenifldNetTable');

% 'Ttenifldxct', 'TenifldXctTable'
% 'Ttenifldnet', 'TenifldNetTable'
% 'Tgseaoutput', 'GSEAResultTable'
if gseaok
    gui.i_exporttable(Tr, true, 'Tgseaoutput', 'GSEAResultTable');
end
end

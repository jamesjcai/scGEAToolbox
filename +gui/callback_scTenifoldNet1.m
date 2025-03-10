function callback_scTenifoldNet1(src, events)
import ten.*
try
    ten.check_tensor_toolbox;
catch ME
    errordlg(ME.message);
    return;
end

extprogname = 'scTenifoldNet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end
if isfolder(wkdir)
    cd(wkdir);
end

%     if exist('sctenifoldnet','file')~=2
%         errordlg('scTenifoldNet is not installed.');
%         disp('To install scTenifoldNet, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
%         disp('addpath(''./scTenifoldNet-master/MATLAB'');');
%         return;
%     end

[FigureHandle, sce, isui] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, 'Construct gene regulatory network (GRN) for all cells or selected cells?', ...
        '', {'All Cells', 'Select Cells...', 'Cancel'}, ...
        'All Cells');
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

    answer = gui.myQuestdlg(FigureHandle, 'This analysis may take several hours. Continue?');
    if ~strcmpi(answer, 'Yes'), return; end
    tmpmat = tempname;
    fw = gui.gui_waitbar;
    try
        disp('>> [A]=ten.sc_pcnetdenoised(sce.X,''savegrn'',false);');
        [A] = ten.sc_pcnetdenoised(sce.X, 'savegrn', false);
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    try
        g = sce.g;
        fprintf('Saving network (A) to %s.mat\n', tmpmat);
        save(tmpmat, 'A', 'g');
    catch ME
        disp(ME);
    end
    % tstr=matlab.lang.makeValidName(string(datetime));
    % save(sprintf('A_%s',tstr),'A','g','-v7.3');
    if ~(ismcc || isdeployed)
        labels = {'Save network to variable named:', ...
            'Save sce.g to variable named:'};
        vars = {'A', 'g'};
        values = {A, sce.g};
        export2wsdlg(labels, vars, values);
    end
end
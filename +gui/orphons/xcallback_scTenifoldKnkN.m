function xcallback_scTenifoldKnkN(src, ~)
import ten.*
%     if exist('sctenifoldnet','file')~=2
%         errordlg('scTenifoldNet is not installed.');
%         disp('To install scTenifoldNet, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
%         disp('addpath(''./scTenifoldNet-master/MATLAB'');');
%         return;
%     end
%     if exist('sctenifoldknk','file')~=2
%         errordlg('scTenifoldKnk is not installed.');
%         disp('To install scTenifoldKnk, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldKnk/archive/master.zip'');');
%         disp('addpath(''./scTenifoldKnk-master/MATLAB'');');
%         return;
%     end

if (ismcc || isdeployed)
    errordlg('This function is not for standalone application.');
            return;
        end

        FigureHandle = src.Parent.Parent;
        sce = guidata(FigureHandle);

        answer = questdlg('Construct network de novo or use existing network in Workspace?', ...
            'Input Network', 'Construct de novo', 'Use existing', 'Construct de novo');
        switch answer
            case 'Use existing'
                a = evalin('base', 'whos');
                b = struct2cell(a);
                if isempty(b)
                    helpdlg('No variable in the WorkSpace.', '');
                    return;
                end
                valididx = ismember(b(4, :), 'double');
                b1 = b(1, :);
                b2 = b(2, :);
                targetsiz = [length(sce.g), length(sce.g)];
                for k = 1:length(valididx)
                    if valididx(k)
                        valididx(k) = isequal(b2{k}, targetsiz);
                    end
                end
                if ~any(valididx)
                    warndlg('No valid network.');
                    return;
                end
                [indx, tf] = listdlg('PromptString', ...
                    {'Select network:'}, ...
                    'liststring', b1(valididx), ...
                    'SelectionMode', 'single');
                if tf ~= 1, return; end
                [A0] = evalin('base', a(indx).name);
                [m, n] = size(A0);
                if m ~= n || n ~= length(sce.g)
                    errordlg('Not a valid network.');
                    return;
                end
            case 'Construct de novo'
                try
                    ten.check_tensor_toolbox;
                catch ME
                    errordlg(ME.message);
                    return;
                end
                A0 = [];
            otherwise
                return;
        end

        answer = questdlg('This analysis may take several hours. Continue?');
        if ~strcmpi(answer, 'Yes'), return; end


        if isempty(A0)
            try
                fw = gui.gui_waitbar;
                fprintf('\nCommand line: [A0]=sc_pcnetdenoised(sce.X);\n');
                fprintf('\nCommand line: [F]=ten.knk3_buildPerturbationLandscape(A0,sce.g);\n');
                [A0] = sc_pcnetdenoised(sce.X);
                [F] = ten.knk3_buildPerturbationLandscape(A0, sce.g);
                gui.gui_waitbar(fw);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
                return;
            end
            isreconstructed = true;
        else
            try
                fw = gui.gui_waitbar;
                fprintf('\nCommand line: [F]=ten.knk3_buildPerturbationLandscape(A,g);\n');
                [F] = ten.knk3_buildPerturbationLandscape(A0, sce.g);

                gui.gui_waitbar(fw);
            catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
                return;
            end
            isreconstructed = false;
        end

        if isreconstructed
            labels = {'Save network to variable named:', ...
                'Save perturbation score to variable named:', ...
                'Save gene list to variable named:'};
            vars = {'A0', 'F', 'g'};
            values = {A0, F, sce.g};
        else
            labels = {'Save perturbation score to variable named:', ...
                'Save gene list to variable named:'};
            vars = {'F', 'g'};
            values = {F, sce.g};
        end
        waitfor(export2wsdlg(labels, vars, values));

    end

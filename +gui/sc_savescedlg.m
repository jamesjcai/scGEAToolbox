function [OKPressed] = sc_savescedlg(sce, parentfig)

if nargin<2, parentfig = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

OKPressed = false;

list = {'SCE Data File (*.mat)...', ...
        'Seurat/Rds File (*.rds)...', ...
        'AnnData/H5ad File (*.h5ad)...', ...
        'Export SCE Data to Workspace...'};

preftagname ='savescedlgindex';
defaultindx = getpref('scgeatoolbox', preftagname, length(list));

       if gui.i_isuifig(parentfig)
            [indx, tf] = gui.myListdlg(parentfig, list, ...
                'Select a destination:', list(defaultindx), false);
        else
            [indx, tf] = listdlg('ListString', list, ...
                'SelectionMode', 'single', ...
                'PromptString', {'Select a destination:'}, ...
                'ListSize', [220, 300], ...
                'Name', 'Export Data', ...
                'InitialValue', defaultindx);
        end
if tf ~= 1, return; end

a = sce.metadata(contains(sce.metadata, "Source:"));
if ~isempty(a), a = strtrim(strrep(a, "Source: ","")); end
if ~isempty(a), a = strrep(a,sprintf("\nOrganism:"),""); end
if ~isempty(a), a = strjoin(a(:)', '_'); end
if ~isempty(a), a = matlab.lang.makeValidName(a); end

setpref('scgeatoolbox', preftagname, indx);
ButtonName = list{indx};
switch ButtonName
        case 'SCE Data File (*.mat)...'
            if ~isempty(a)
                [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as');
            end
            if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(parentfig);
                % Disk full, read-only folder or a locked file: the bar stayed
                % up and the menu callback errored.
                try
                    save(filename, 'sce', '-v7.3');
                catch ME
                    gui.myWaitbar(parentfig, fw, true);
                    gui.myErrordlg(parentfig, ME.message, ME.identifier);
                    return;
                end
                gui.myWaitbar(parentfig, fw);
                OKPressed = true;
            end
        case 'Seurat/Rds File (*.rds)...'
            % One question, and only when R is missing. It used to ask
            % "requires R. Continue?" every time, then after the file
            % picker ask again whether to set R up.
            if ~ispref('scgeatoolbox', 'rexecutablepath')
                if strcmp(gui.myQuestdlg(parentfig, ['Saving as .rds needs R, ' ...
                        'which is not set up. Set it up now?']), 'Yes')
                    gui.i_setrenv(parentfig);
                end
                if ~ispref('scgeatoolbox', 'rexecutablepath'), return; end
            end
            if ~isempty(a)
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as');
            end
            if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(parentfig);
                % A failed save used to be reported as a success: the
                % status was ignored and no file existed.
                try
                    status = sc_sce2rds(sce, filename);
                catch ME
                    gui.myWaitbar(parentfig, fw, true);
                    gui.myErrordlg(parentfig, ME.message, ME.identifier);
                    return;
                end
                if ~status || ~isfile(filename)
                    gui.myWaitbar(parentfig, fw, true);
                    gui.myErrordlg(parentfig, sprintf(['R did not write %s. ' ...
                        'The R output in the Command Window should say why.'], filename));
                    return;
                end
                gui.myWaitbar(parentfig, fw);
                fprintf("\nTo read file, in R:\n");
                fprintf("library(Seurat)\n");
                fprintf("A<-readRDS(""%s"")\n", file);
                OKPressed = true;
            end
        case 'AnnData/H5ad File (*.h5ad)...'
            % No Python question: the .h5ad writer is native MATLAB now.
            if ~isempty(a)
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
            end
            if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                try
                    ok = sc_sce2h5ad(sce, filename);
                catch ME
                    % It rethrows after removing the partial file.
                    gui.myErrordlg(parentfig, ME.message, ME.identifier);
                    return;
                end
                if ok
                    fprintf("\nTo read file, in Python:\n");
                    fprintf("adata = anndata.read_h5ad(""%s"")\n", file);
                    OKPressed = true;
                end
            end

        case 'Export SCE Data to Workspace...'
            labels = {'Save SCE to variable named:', ...
                'Save SCE.X to variable named:', ...
                'Save SCE.g to variable named:', ...
                'Save SCE.S to variable named:'};
            vars = {'sce', 'X', 'g', 's'};
            values = {copy(sce), sce.X, sce.g, sce.s};

            if gui.i_isuifig(parentfig)
                [~, OKPressed] = gui.myExport2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', ...
                    [true, false, false, false], parentfig);
            else
                OKPressed = gui.i_export2wsdlg(parentfig, labels, vars, values, ...
                    'Save Data to Workspace', ...
                    logical([1, 0, 0, 0]));
            end
        otherwise
            return;
    end
end

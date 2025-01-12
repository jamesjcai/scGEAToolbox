function [OKPressed] = sc_savescedlg(sce)

    OKPressed = false;

    list = {'SCE Data File (*.mat)...', ...
        'Seurat/Rds File (*.rds)...', ...
        'AnnData/H5ad File (*.h5ad)...', ...
        'Export SCE Data to Workspace...'};
 
    [indx, tf] = listdlg('ListString', list, ...
        'SelectionMode', 'single', ...
        'PromptString', {'Select a destination:'}, ...
        'ListSize', [220, 300], ...
        'Name', 'Export Data', ...
        'InitialValue', length(list));
    if tf ~= 1, return; end

            a = sce.metadata(contains(sce.metadata, "Source:"));
            if ~isempty(a), a = strtrim(strrep(a, "Source: ","")); end
            if ~isempty(a), a = sprintf("%s_", a(:)); end
            if ~isempty(a), a = extractBefore(a, strlength(a)); end
            if ~isempty(a), a = matlab.lang.makeValidName(a); end
    ButtonName = list{indx};
    switch ButtonName
        case 'SCE Data File (*.mat)...'
            if ~isempty(a)
                [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as');
            end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.gui_waitbar;
                save(filename, 'sce', '-v7.3');
                gui.gui_waitbar(fw);
                OKPressed = true;
            end
        % case 'TXT/TSV/CSV File (*.txt)...'
        %     warndlg('Function is under development.');
        %     return;
        %     if ~isempty(a)
        %         [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as', a);
        %     else
        %         [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
        %     end
        %     if isequal(file, 0) || isequal(path, 0)
        %         return;
        %     else
        %         filename = fullfile(path, file);
        %         fw = gui.gui_waitbar;
        % 
        %         gui.gui_waitbar(fw);
        %         OKPressed = true;
        %     end            
        case 'Seurat/Rds File (*.rds)...'
            answer = questdlg('This function requires R. Continue?','');
            if ~strcmp(answer,'Yes'), return; end            
            if ~isempty(a)
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as');
            end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.gui_waitbar;
                sc_sce2rds(sce, filename);
                gui.gui_waitbar(fw);
                fprintf("\nTo read file, in R:\n");
                fprintf("library(Seurat)\n");
                fprintf("A<-readRDS(""%s"")\n", file);
                OKPressed = true;
            end
        case 'AnnData/H5ad File (*.h5ad)...'
            answer = questdlg('This function requires Python. Continue?','');
            if ~strcmp(answer,'Yes'), return; end
            if ~isempty(a)
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
            end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                if sc_sce2h5ad(sce, filename)
                    fprintf("\nTo read file, in Python:\n");
                    fprintf("adata = anndata.read(""%s"")\n", file);
                    OKPressed = true;
                end
            end
            
        case 'Export SCE Data to Workspace...'
            labels = {'Save SCE to variable named:', ...
                'Save SCE.X to variable named:', ...
                'Save SCE.g to variable named:', ...
                'Save SCE.S to variable named:'};
            vars = {'sce', 'X', 'g', 's'};
            values = {sce, sce.X, sce.g, sce.s};
            [~, OKPressed] = export2wsdlg(labels, vars, values, ...
                'Save Data to Workspace', ...
                logical([1, 0, 0, 0]), {@smhelp});
        otherwise
            return;
    end
end

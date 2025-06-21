function [OKPressed] = callback_SaveX(app, events)
    OKPressed = false;

    list = {'SCE Data File (*.mat)...', ...
        'Seurat/Rds File (*.rds)...', ...
        'AnnData/H5ad File (*.h5ad)...', ...
        'Export SCE Data to Workspace...'};

    preftagname ='savescedlgindex';
    defaultindx = getpref('scgeatoolbox', preftagname, length(list));

    [indx, tf] = gui.myListdlg(app.UIFigure, list, ...
        'Select a destination:', list(defaultindx));
    if tf ~= 1, return; end
            a = app.sce.metadata(contains(app.sce.metadata, "Source:"));
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
            if isvalid(app.UIFigure) && isa(app.UIFigure, 'matlab.ui.Figure'), figure(app.UIFigure); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(app.UIFigure);
                sce = app.sce;                
                save(filename, 'sce', '-v7.3');
                gui.myWaitbar(app.UIFigure, fw);
                OKPressed = true;
            end
        case 'Seurat/Rds File (*.rds)...'
            answer = gui.myQuestdlg(app.UIFigure, 'This function requires R. Continue?','');
            if ~strcmp(answer,'Yes'), return; end            
            if ~isempty(a)
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as');
            end
            if isvalid(app.UIFigure) && isa(app.UIFigure, 'matlab.ui.Figure'), figure(app.UIFigure); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(app.UIFigure);
                sc_sce2rds(app.sce, filename);
                gui.myWaitbar(app.UIFigure, fw);
                fprintf("\nTo read file, in R:\n");
                fprintf("library(Seurat)\n");
                fprintf("A<-readRDS(""%s"")\n", file);
                OKPressed = true;
            end
        case 'AnnData/H5ad File (*.h5ad)...'
            answer = gui.myQuestdlg(app.UIFigure, 'This function requires Python. Continue?','');
            if ~strcmp(answer,'Yes'), return; end
            if ~isempty(a)
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
            end
            if isvalid(app.UIFigure) && isa(app.UIFigure, 'matlab.ui.Figure'), figure(app.UIFigure); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                if sc_sce2h5ad(app.sce, filename)
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
            values = {app.sce, app.sce.X, app.sce.g, app.sce.s};
            if gui.i_isuifig(app.UIFigure)
                [~, OKPressed] = gui.myExport2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', app.UIFigure);
            else
                [~, OKPressed] = export2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', ...
                    logical([1, 0, 0, 0]), {@smhelp});
            end
        otherwise
            return;
    end
end
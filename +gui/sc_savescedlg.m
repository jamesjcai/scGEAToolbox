function [OKPressed] = sc_savescedlg(sce, parentfig)

if nargin<2, parentfig = []; end

    OKPressed = false;

    list = {'SCE Data File (*.mat)...', ...
        'Seurat/Rds File (*.rds)...', ...
        'AnnData/H5ad File (*.h5ad)...', ...
        'Export SCE Data to Workspace...'};

    preftagname ='savescedlgindex';
    defaultindx = getpref('scgeatoolbox', preftagname, length(list));

       if gui.i_isuifig(parentfig)
            [indx, tf] = gui.myListdlg(parentfig, list, ...
                'Select a destination:', list(defaultindx));
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
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(parentfig);
                save(filename, 'sce', '-v7.3');
                gui.myWaitbar(parentfig, fw);
                OKPressed = true;
            end
        case 'Seurat/Rds File (*.rds)...'
            answer = gui.myQuestdlg(parentfig, 'This function requires R. Continue?','');
            if ~strcmp(answer,'Yes'), return; end            
            if ~isempty(a)
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.rds'; '*.*'}, 'Save as');
            end
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(parentfig);
                sc_sce2rds(sce, filename);
                gui.myWaitbar(parentfig, fw);
                fprintf("\nTo read file, in R:\n");
                fprintf("library(Seurat)\n");
                fprintf("A<-readRDS(""%s"")\n", file);
                OKPressed = true;
            end
        case 'AnnData/H5ad File (*.h5ad)...'
            answer = gui.myQuestdlg(parentfig, 'This function requires Python. Continue?','');
            if ~strcmp(answer,'Yes'), return; end
            if ~isempty(a)
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as', a);
            else
                [file, path] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
            end
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
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
            if gui.i_isuifig(parentfig)
                [~, OKPressed] = gui.myExport2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', ...
                    [true, false, false, false], parentfig);
            else
                [~, OKPressed] = export2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', ...
                    logical([1, 0, 0, 0]), {@smhelp});
            end
        otherwise
            return;
    end
end

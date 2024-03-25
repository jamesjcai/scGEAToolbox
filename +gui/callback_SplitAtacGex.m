function [sceg, scep] = callback_SplitAtacGex(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

if ~any(contains(sce.g, ':'))
    warndlg('Not a multiome ATAC+GEX matrix.');
    return;
end

ispeak = contains(sce.g, ':');
isgene = ~ispeak;

sceg = sce;
sceg.g = sceg.g(isgene);
sceg.X = sceg.X(isgene, :);

scep = sce;
scep.g = scep.g(~isgene);
scep.X = scep.X(~isgene, :);


labels = {'Save SCEG to variable named:', ...
    'Save SCEP to variable named:'};
vars = {'sceg', 'scep'};
values = {sceg, scep};
[~, ~] = export2wsdlg(labels, vars, values, ...
    'Save Data to Workspace');

% --------------------
pause(2);
fx = scgeatool(sceg);
fx.Position(3:4) = 0.8 * fx.Position(3:4);
movegui(fx, 'center');
fx.Position(1) = fx.Position(1) - 250;
fx = fx.CurrentAxes;
fy.Subtitle.String = '[genes x cells]';
%waitfor(helpdlg(sprintf('%s Cells extracted.', ...
%    sprintf('%s+',tg)),''));
answer = questdlg('scRNAseq data extracted. Continue?', '');
if ~strcmp(answer, 'Yes')
    return;
end
fy = scgeatool(scep);
fy.Position(3:4) = 0.8 * fy.Position(3:4);
movegui(fy, 'center');
fy.Position(1) = fy.Position(1) + 250;
fy = fy.CurrentAxes;
fy.Subtitle.String = '[peaks x cells]';


% pause(4);
% scgeatool(sceg);
%
% pause(4);
% scgeatool(scep);
% ---------------------

%{
%if ~(ismcc || isdeployed)
answer = questdlg('Export & save data to:','',...
    'Workspace','MAT file','Seurat/RDS file','Workspace');
%else
%    answer = questdlg('Export & save data to:','',...
%        'MAT file','Seurat/RDS file','MAT file');
%end

switch answer
    case 'Workspace'
        labels = {'Save SCE to variable named:',...
            'Save SCE.X to variable named:',...
            'Save SCE.g to variable named:',...
            'Save SCE.S to variable named:'};
        vars = {'sce','X','g','s'};
        values = {sce,sce.X,sce.g,sce.s};
        [~,OKPressed]=export2wsdlg(labels,vars,values,...
            'Save Data to Workspace',...
            logical([1 0 0 0]),{@smhelp});
    case 'MAT file'
        [file, path] = uiputfile({'*.mat';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
            return;
        else
            filename=fullfile(path,file);
            fw=gui.gui_waitbar;
            save(filename,'sce','-v7.3');
            gui.gui_waitbar(fw);
            OKPressed=true;
        end
    case 'Seurat/RDS file'
        [file, path] = uiputfile({'*.rds';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
            return;
        else
            filename=fullfile(path,file);
            fw=gui.gui_waitbar;
            sc_sce2rds(sce,filename);
            gui.gui_waitbar(fw);
            disp("A<-readRDS(""input.rds"")");
            OKPressed=true;
        end
    otherwise
        return;
end
    function smhelp
        helpdlg({'Select one or both check boxes.',...
            'Change the variable names, if desired,',...
                'and then click OK.'});
        end
end
    %}
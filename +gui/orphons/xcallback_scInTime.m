function xcallback_scInTime(src, ~)
import ten.*
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

ptime = sce.list_cell_attributes{2};

answer = questdlg('This analysis may take several hours. Continue?');
if ~strcmpi(answer, 'Yes'), return; end
tmpmat = tempname;
fw = gui.gui_waitbar;
try
    disp('>> [T]=scInTime(sce.X,sce.g,ptime);');
    [s] = scInTime(sce.X, sce.g, ptime);
catch ME
    gui.gui_waitbar(fw);
    errordlg(ME.message);
    return;
end
gui.gui_waitbar(fw);

%     if ~(ismcc || isdeployed)
%         labels = {'Save network to variable named:',...
%             'Save sce.g to variable named:'};
%         vars = {'A','g'};
%         values = {A,sce.g};
%         export2wsdlg(labels,vars,values);
%     end

end
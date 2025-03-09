function [wrkdir] = gui_setprgmwkdir(extprogname, preftagname, FigureHandle)
if nargin<3, FigureHandle = []; end
wrkdir = '';
%extprogname = 'R_monocle3';
%preftagname = 'externalwrkpath';
if ~gui.i_setwrkdir(preftagname, FigureHandle), return; end
s = getpref('scgeatoolbox', preftagname, []);
if isempty(s)
    error('Working path has not been set up.');
end
s1 = sprintf('%s_workingfolder', extprogname);
wrkdir = fullfile(s, s1);

if ~exist(wrkdir,"dir")
    mkdir(wrkdir);
else
    answer = gui.myQuestdlg(FigureHandle, sprintf('%s existing. Overwrite?', wrkdir));
    if ~strcmp(answer,'Yes')
        wrkdir = '';
        return;
    else
        if ~strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Existing files in the working folder will be overwritten or deleted. Continue?'))
            wrkdir = '';
            return;          
        end
    end
end
   
fprintf('CURRENTWDIR = "%s"\n', wrkdir);

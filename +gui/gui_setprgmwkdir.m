function [wrkdir] = gui_setprgmwkdir(extprogname, preftagname)

wrkdir = '';
%extprogname = 'R_monocle3';
%preftagname = 'externalwrkpath';
if ~gui.i_setwrkdir(preftagname), return; end
s = getpref('scgeatoolbox', preftagname, []);
if isempty(s)
    error('Working path has not been set up.');
end
s1 = sprintf('%s_workingfolder', extprogname);
wrkdir = fullfile(s, s1);

if ~exist(wrkdir,"dir")
    mkdir(wrkdir);
else
    answer = questdlg(sprintf('%s existing. Overwrite?', wrkdir));
    if ~strcmp(answer,'Yes')
        wrkdir = '';
        return;
    else
        if ~strcmp('Yes', questdlg('Existing files in the working folder will be overwritten or deleted. Continue?'))
            wrkdir = '';
            return;          
        end
    end
end
   
fprintf('CURRENTWDIR = "%s"\n', wrkdir);

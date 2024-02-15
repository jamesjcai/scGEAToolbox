function [wkdir] = gui_setprgmwkdir(extprogname, preftagname)

wkdir = '';
%extprogname = 'R_monocle3';
%preftagname = 'externalwrkpath';
if ~gui.i_setwrkdir(preftagname), return; end
s = getpref('scgeatoolbox', preftagname,[]);
if isempty(s)
    error('Working path has not been set up.');
end
s1 = sprintf('%s_workingfolder', extprogname);
wkdir = fullfile(s, s1);

if ~exist(wkdir,"dir")
    mkdir(wkdir);
else
    answer = questdlg(sprintf('%s existing. Overwrite?', ...
        wkdir));
    if ~strcmp(answer,'Yes')
        wkdir = '';
        return;
    end
end
   
fprintf('CURRENTWDIR = "%s"\n', wkdir);
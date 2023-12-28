function [wkdir] = gui_setprgmwkdir(extprogname, preftagname)

wkdir = '';
%extprogname = 'R_monocle3';
%preftagname = 'externalwrkpath';
if ~gui.i_setwrkdir(preftagname), return; end
s = getpref('scgeatoolbox', preftagname);
s1 = sprintf('%s_workingfolder', extprogname);
wkdir = fullfile(s, s1);

if ~exist(wkdir,"dir")
    mkdir(wkdir);
else
    answer = questdlg('Directory existing. Overwrite?');
    if ~strcmp(answer,'Yes'), return; end
end
   
fprintf('CURRENTWDIR = "%s"\n', wkdir);
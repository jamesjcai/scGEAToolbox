function [ok, msg, codepth] = commoncheck_R(rscriptdir, externalfolder)

if nargin < 2, externalfolder = 'external'; end
ok = false;
msg = [];

if ~ispref('scgeatoolbox', 'rexecutablepath')
    answer = questdlg('Select R Interpreter?');
    if strcmp(answer, 'Yes'), gui.i_setrenv; end
    return;
else
    Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
    if isempty(Rpath)
        error('R environment has not been set up.');
    end    
end

%if isempty(pkg.FindRpath)
%   msg=('Rscript.exe is not found');
%   return;
%end


folder = fileparts(mfilename('fullpath'));
a = strfind(folder, filesep);
folder = extractBefore(folder, a(end)+1);
codepth = fullfile(folder, externalfolder, rscriptdir);
cd(codepth);
% fprintf('CURRENTWDIR = "%s"\n', wrkpth);
fprintf('SORCCODEDIR = "%s"\n', codepth);
[~, cmdout] = pkg.RunRcode('require.R', Rpath);
if strfind(cmdout, 'there is no package') > 0
    msg = cmdout;
    return;
end
ok = true;
% disp('commoncheck_R is Done.');
end

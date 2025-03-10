function [ok, msg] = i_commoncheck_R(extprogname, FigureHandle)

if nargin<2, FigureHandle = []; end
if nargin < 1, extprogname = 'r_enrichR'; end

ok = false;
if ~ispref('scgeatoolbox', 'rexecutablepath')
    answer = gui.myQuestdlg(FigureHandle, 'Select R Interpreter?');
    if strcmp(answer, 'Yes'), gui.i_setrenv; end
    return;
else
    Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
    if isempty(Rpath)
        warning('R environment has not been set up.');
        return;
    end    
end
oldpath = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', '+run', 'external', extprogname);
cd(codepth);

[~, cmdout] = pkg.RunRcode('require.R', Rpath);
cd(oldpath);
if strfind(cmdout, 'there is no package') > 0
    msg = cmdout;
    return;
end
    ok = true;
end
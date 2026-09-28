function [ok, msg, codepth] = commoncheck_R(rscriptdir, externalfolder, FigureHandle)

if nargin < 3, FigureHandle = []; end
if nargin < 2, externalfolder = 'external'; end
ok = false;
% Every output assigned up front. CODEPTH was left unset on the early
% returns, so all 16 callers, which ask for three outputs, threw "Output
% argument not assigned"; and an empty MSG made their `error(msg)` a no-op.
codepth = '';
msg = 'R has not been set up. Set it with Setup > Set Up R Environment, then try again.';

if ~ispref('scgeatoolbox', 'rexecutablepath')
    answer = gui.myQuestdlg(FigureHandle, 'Select R Interpreter?');
    if strcmp(answer, 'Yes'), gui.i_setrenv(FigureHandle); end
    % Carry on if the user has just set it up.
    if ~ispref('scgeatoolbox', 'rexecutablepath'), return; end
end
Rpath = getpref('scgeatoolbox', 'rexecutablepath', []);
if isempty(Rpath), return; end
msg = [];


folder = fileparts(mfilename('fullpath'));
a = strfind(folder, filesep);
folder = extractBefore(folder, a(end)+1);


codepth = fullfile(folder, externalfolder, rscriptdir);
if ~exist(codepth,"dir")
    codepth = fullfile(folder, '..', externalfolder, rscriptdir);
end
if ~exist(codepth,"dir")
    error('CODEPTH is undefined.');
end


codepth = pkg.i_normalizepath(codepth);

fprintf('R___CODEDIR = "%s"\n', codepth);
[~, cmdout] = pkg.i_runrcode(fullfile(codepth, 'require.R'), Rpath);
if ~isempty(strfind(cmdout, 'there is no package'))
    msg = cmdout;
    return;
end
ok = true;
end

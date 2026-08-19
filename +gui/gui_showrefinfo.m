function [y, txt, T] = gui_showrefinfo(reftarget, parentfig)
% see also: gui.gui_uishowrefinfo

if nargin<2, parentfig = []; end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end

y=false;
txt = [];
pw1 = fileparts(mfilename('fullpath'));
fname = fullfile(pw1, '..','assets','Misc','refinfo.txt');
fid=fopen(fname,'r');
T=textscan(fid,'%s%s','Delimiter','\t');
fclose(fid);
reftag=string(T{:,1});

idx=find(reftarget==reftag);
if ~isempty(idx)
    txt=T{:,2}{idx};
end

if isempty(txt) , return; end

% refinfo.txt stores line breaks as the literal characters "\n" (and tabs as
% "\t"). Convert them to real newline/tab characters so the dialog renders on
% multiple lines instead of showing "\n". Use literal replacement rather than
% compose/sprintf so any "%" in the help text is left untouched.
txt = replace(string(txt), "\t", sprintf("\t"));
txt = replace(txt, "\n", newline);

% Drop a trailing ";" (with any trailing spaces) at the end of each line.
% These were list separators that are redundant now that each item renders
% on its own line.
txt = regexprep(txt, ';[ \t]*$', '', 'lineanchors');

if nargout == 0
    % Informational display (e.g. Help > Shortcuts User Guide). A single
    % dismiss button is the right affordance here, not Continue/Cancel.
    gui.myHelpdlg(parentfig, txt, reftarget);
    y = true;
else
    % Used as a confirmation gate by callers, e.g.:
    %   if ~gui.gui_showrefinfo(target, fig), return; end
    answer = gui.myQuestdlg(parentfig, ...
        txt, reftarget, {'Continue', 'Cancel'}, 'Continue');
    if strcmp(answer, 'Continue')
        y = true;
    end
end
end

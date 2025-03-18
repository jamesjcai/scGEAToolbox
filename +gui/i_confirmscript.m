function [t] = i_confirmscript(qtxt, stxt, langtag, parentfig)

if nargin<4, parentfig = []; end

t = false;
if nargin < 1, qtxt = 'Run pseudotime analysis (Monocle)?'; end
if nargin < 2, stxt = 'R_monocle'; end
if nargin < 3, langtag = 'R'; end

switch lower(langtag)
    case 'r'
        scriptfile = 'script.R';
    case 'python'
        scriptfile = 'script.py';
    otherwise
        error('Unknown language tag.');
end

answer = gui.myQuestdlg(parentfig, qtxt, '', ...
    {'Yes', 'Review Script', 'Cancel'}, 'Yes');
switch answer
    case 'Cancel'
        return;
    case 'Yes'
        t = true;
    case 'Review Script'
        folder = fileparts(mfilename('fullpath'));
        scriptfile = fullfile(folder, '..', '+run', 'external', ...
            stxt, scriptfile);
        ts = fileread(scriptfile);
        %LF=char(10);
        CR = char(13); %  carriage return character equivalent to char(13) or sprintf('\r').
        ts = strrep(ts, [CR, newline], newline);

        if gui.i_isuifig(parentfig)
            a = gui.myInputdlg({'Review script and press OK to run it'}, '', {ts}, parentfig);
        else
            a = inputdlg('Review script and press OK to run it', '', [15, 80], {ts});
        end

        if isempty(a)
            return;
        else
            answer2 = gui.myQuestdlg(parentfig, "Run script?", '');
            if isempty(answer2), return; end
            switch answer2
                case 'Yes'
                    t = true;
                otherwise
                    return;
            end
        end
    otherwise
        return;
end

end

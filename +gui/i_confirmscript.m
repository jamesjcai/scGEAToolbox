function [t]=i_confirmscript(qtxt,stxt,langtag)

t=false;
if nargin<1, qtxt='Run pseudotime analysis (Monocle)?'; end
if nargin<2, stxt='R_monocle'; end
if nargin<3, langtag='R'; end

switch lower(langtag)
    case 'r'
        scriptfile='script.R';
    case 'python'
        scriptfile='script.py';
    otherwise
        error('Unknown language tag.');
end

        answer = questdlg(qtxt,'', ...
            'Yes','Review Script','Cancel','Yes');
        switch answer
            case 'Cancel'
                return;
            case 'Yes'
                t=true;
             case 'Review Script'
                folder=fileparts(mfilename('fullpath'));
                scriptfile=fullfile(folder,'..','+run','external', ...
                    stxt,scriptfile);
                ts=fileread(scriptfile);
                %LF=char(10);
                CR=char(13);  %  carriage return character equivalent to char(13) or sprintf('\r').
                ts=strrep(ts,[CR newline],newline);
                a=inputdlg('Review script and press OK to run it', ...
                    '',[10 90],{ts});
                if isempty(a)
                    return; 
                else
                    t=true;
                end
            otherwise
                return;
        end

end

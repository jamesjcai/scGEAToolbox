function [status, cmdout] = RunRcode(RscriptFileName, Rexecpath)
% This function calls R to run R script (.r file) under Matlab, and returns 'runR.log' in the same folder with the input R script file.
% This code only works in Windows environments. It might work in Mac by modifying 'FindRpath'.
% 'RscriptFileName' : path + filename of the R script to be run.
% 'Rpath' (optional) : the path for the installed 'R.exe'.  e.g.  Rpath = 'C:\Program Files\R\R-3.1.1\bin';
% If 'Rpath' is not provided, 'FindRpath' will be executed to automatically find the path for R in Windows Program Files folder.
% Example:
% >> Rpath = 'C:\Program Files\R\R-3.1.1\bin';
% >> RscriptFileName = 'D:\test\lmm.R';
% >> RunRcode(RscriptFileName, Rpath);
% Update:
% Ver. 1.4  Dec-14-2017  support parallel computing (run several R codes simultaneously)
% Weirong Chen   March-8-2015

status = 99;
cmdout = [];

if ~exist(RscriptFileName, 'file')
    disp('RscriptFile not found.');
    return;
end
narginchk(2, 2)

if ispc
    Rexec = fullfile(Rexecpath, 'Rscript.exe');
    if ~exist(Rexec, 'file')
        error([Rexec, ' not existing.'])
    end
else
    Rexec = fullfile(Rexecpath, 'Rscript');
    if ~exist(Rexec, 'file')
        Rexec = fullfile(Rexecpath, 'R');
        if ~exist(Rexec, 'file')
            error([Rexec, ' not existing.'])
        end
    end
end


% if nargin<2 || isempty(Rpath), Rpath=pkg.FindRpath; end
% if isempty(Rpath)
%     if ~ispref('scgeatoolbox','rexecutablepath')
%         return;
%     else
%         Rpath=getpref('scgeatoolbox','rexecutablepath',[]);
% if isempty(Rpath)
%     error('R environment has not been set up.');
% end
%     end
% end

% sep=filesep;
% [p,f,~]=fileparts(RscriptFileName);
% if isempty(p), p = pwd; end
% logFName=[p sep f '.R.log'];

% if iscell(Rpath)
%     Rpath=Rpath{end};
% end

if ispc
    % commandline=['"' Rpath sep 'R.exe" CMD BATCH "' RscriptFileName '" "' logFName '"'];
    % commandline=['"' Rpath sep 'Rscript.exe" "' RscriptFileName '"'];
    commandline = ['"', Rexec, '" "', RscriptFileName, '"'];
    fprintf('COMMANDLINE = %s\n', commandline);
    [status, cmdout] = system(commandline, '-echo');
else
    %commandline=[Rpath ' ' RscriptFileName];
    %fprintf('COMMANDLINE = %s\n',[Rpath ' ' RscriptFileName]);
    %commandline=[Rpath sep 'Rscript ' RscriptFileName];
    if strcmp(Rexec(end-6:end), 'Rscript')
        % commandline=[Rexec ' ' RscriptFileName];
        commandline = sprintf("%s %s", Rexec, RscriptFileName);
        fprintf('COMMANDLINE = %s\n', commandline);
        [status, cmdout] = system(commandline);
    else
        %commandline=[Rexec ' CMD BATCH ' RscriptFileName];
        commandline = sprintf("%s CMD BATCH %s", Rexec, RscriptFileName);
        fprintf('COMMANDLINE = %s\n', commandline);
        [status, cmdout] = system(commandline);
    end
end

end %RunRcode

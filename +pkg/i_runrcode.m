function [status, cmdout] = i_runrcode(RscriptFileName, Rexecpath)
% This function calls R to run R script (.r file) under Matlab, and returns 'runR.log' in the same folder with the input R script file.
% This code only works in Windows environments. It might work in Mac by modifying 'i_findrpath'.
% 'RscriptFileName' : path + filename of the R script to be run.
% 'Rpath' (optional) : the path for the installed 'R.exe'.  e.g.  Rpath = 'C:\Program Files\R\R-3.1.1\bin';
% If 'Rpath' is not provided, 'i_findrpath' will be executed to automatically find the path for R in Windows Program Files folder.
% Example:
% >> Rpath = 'C:\Program Files\R\R-3.1.1\bin';
% >> RscriptFileName = 'D:\test\lmm.R';
% >> i_runrcode(RscriptFileName, Rpath);
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

Rexec = i_resolve_rexec(Rexecpath);
if ~exist(Rexec, 'file')
    error([Rexec, ' not existing.'])
end

[~, rexecName, rexecExt] = fileparts(Rexec);
rexecTag = lower([rexecName, rexecExt]);
if strcmp(rexecTag, 'r') || strcmp(rexecTag, 'r.exe')
    commandline = sprintf('"%s" CMD BATCH "%s"', Rexec, RscriptFileName);
else
    commandline = sprintf('"%s" "%s"', Rexec, RscriptFileName);
end

fprintf('COMMANDLINE = %s\n', commandline);
if ispc
    [status, cmdout] = system(commandline, '-echo');
else
    [status, cmdout] = system(commandline);
end

end % i_runrcode

function Rexec = i_resolve_rexec(Rexecpath)
if isfolder(Rexecpath)
    if ispc
        Rexec = fullfile(Rexecpath, 'Rscript.exe');
    else
        Rexec = fullfile(Rexecpath, 'Rscript');
        if ~exist(Rexec, 'file')
            Rexec = fullfile(Rexecpath, 'R');
        end
    end
else
    Rexec = char(Rexecpath);
end
end

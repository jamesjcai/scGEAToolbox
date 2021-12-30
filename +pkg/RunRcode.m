function [status,cmdout]=RunRcode(RscriptFileName,Rpath)
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

status=99; 
cmdout=[];
if ~exist(RscriptFileName,'file'), return; end
if nargin<2 || isempty(Rpath), Rpath=pkg.FindRpath; end
if isempty(Rpath), return; end

sep=filesep;
% [p,f,~]=fileparts(RscriptFileName);
% if isempty(p), p = pwd; end
% logFName=[p sep f '.R.log'];

if iscell(Rpath)
    Rpath=Rpath{end};    
end

if ispc
    % commandline=['"' Rpath sep 'R.exe" CMD BATCH "' RscriptFileName '" "' logFName '"'];
    commandline=['"' Rpath sep 'Rscript.exe" "' RscriptFileName '"'];
    fprintf('COMMANDLINE = %s\n',commandline);
elseif isunix
    commandline=[Rpath ' ' RscriptFileName];
    fprintf('COMMANDLINE = %s\n',[Rpath ' ' RscriptFileName]);
end
[status,cmdout]=system(commandline,'-echo');

end %RunRcode


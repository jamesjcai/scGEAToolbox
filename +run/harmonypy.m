function [s]=harmonypy(s,batchid)

if nargin<2, error('[s]=run.harmonypy(s,batchid)'); end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'thirdparty','harmony');
cd(wrkpth);

if exist('output.csv','file'), delete('output.csv'); end
writematrix(s,'input1.csv');
batchidx=matlab.lang.makeValidName(string(batchid));
writetable(table(batchidx),'input2.csv','QuoteStrings',true);
% fid=fopen('input2.csv','w');
% fprintf(fid,'"batchid"\n');
% for k=1:length(batchid)
%     fprintf(fid,'"ID%d"\n',batchid(k));
% end
% fclose(fid);
% pyenv('Version','d:\\Miniconda3\\envs\\harmonypy\\python.exe')
x=pyenv;
cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
disp(cmdlinestr)
system(cmdlinestr);
if exist('output.csv','file')
    s=readmatrix('output.csv');
else    
    % s=[];
end
if exist('input1.csv','file'), delete('input1.csv'); end
if exist('input2.csv','file'), delete('input2.csv'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end


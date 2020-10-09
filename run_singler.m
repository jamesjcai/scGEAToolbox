function [c]=run_singler(X,genelist,species)
if nargin<3, species="human"; end
if nargin<2, genelist=(1:size(X,1))'; end

if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_SingleR');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('output.csv','file')
    delete('output.csv');
end
sc_writefile('input.txt',X,upper(genelist));
switch lower(species)
    case "human"
        disp("human")
        RunRcode('script.R');
    case "mouse"
        RunRcode('script_mouse.R');
        disp("mouse")
end

if exist('output.csv','file')
    T=readtable('output.csv','ReadVariableNames',false);
    c=string(T.Var1);
else
    c=[];
end
%if exist('input.txt','file'), delete('input.txt'); end
%if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);

function [c]=SingleR(X,genelist,species)

isdebug=false;
if nargin<3, species="human"; end
if nargin<2, genelist=(1:size(X,1))'; end

oldpth=pwd();
[isok,msg]=commoncheck_R('R_SingleR');
if ~isok, error(msg); c=[]; return; end

if ~isdebug
if exist('./output.csv','file'), delete('./output.csv'); end
if exist('./input.txt','file'), delete('./input.txt'); end
end

sc_writefile('input.txt',X,upper(genelist));
Rpath=getpref('scgeatoolbox','rexecutablepath');
switch lower(species)
    case "human"
        disp("human")
        pkg.RunRcode('script_human.R',Rpath);
    case "mouse"
        pkg.RunRcode('script_mouse.R',Rpath);
        disp("mouse")
end

if exist('output.csv','file')
    T=readtable('output.csv','ReadVariableNames',false);
    c=string(T.Var1);
else
    c=[];
end
if ~isdebug
if exist('./output.csv','file'), delete('./output.csv'); end
if exist('./input.txt','file'), delete('./input.txt'); end
end
cd(oldpth);
end
function [X2]=SeuratSctransform(X,genelist)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SeuratSctransform');

if ~isok, error(msg); X2=[]; return; end

if exist('output.mat.tmp','file'), delete('output.mat.tmp'); end
%if exist('tsneoutput.csv','file'), delete('tsneoutput.csv'); end
%if exist('umapoutput.csv','file'), delete('umapoutput.csv'); end
%if exist('activeidentoutput.csv','file'), delete('activeidentoutput.csv'); end
%sc_writefile('input.txt',sce.X,sce.g);
if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
    if exist('./input.txt','file'), delete('./input.txt'); end
end

if ~iscellstr(genelist) && isstring(genelist)
    genelist=cellstr(genelist);
end
lastwarn('')
save('input.mat','X','genelist','-v6');
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)    
    disp(warnId)
	if exist('./input.mat','file'), delete('./input.mat'); end
    disp('Writing data into input.txt...')
    sc_writefile('input.txt',X,genelist);
end

pkg.RunRcode('script.R');

if exist('output.mat','file')
    load output.mat X2
end
if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
    if exist('./input.txt','file'), delete('./input.txt'); end
end
cd(oldpth);
end
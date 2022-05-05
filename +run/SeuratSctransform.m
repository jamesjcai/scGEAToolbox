function [X2]=SeuratSctransform(X,genelist)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SeuratSctransform');

if ~isok, error(msg); X2=[]; return; end

pkg.i_deletefiles('output.mat.tmp');
if ~isdebug
    pkg.i_deletefiles('input.mat','output.mat', ...
                      'input.txt','output.txt');
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
elseif exist('output.txt','file')
    X2=readmatrix('output.txt');
end
if ~isdebug
    pkg.i_deletefiles('input.mat','output.mat', ...
                      'input.txt','output.txt');
end
cd(oldpth);
end

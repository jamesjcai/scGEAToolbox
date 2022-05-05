function [X]=SeuratSctransformnew(X,genelist)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SeuratSctransform');
if ~isok, error(msg); return; end
tmpfilelist={'input.mat','output.mat', ...
                      'input.txt','output.txt','g.txt'};

pkg.i_deletefiles('output.mat.tmp');
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
% if ~iscellstr(genelist) && isstring(genelist)
%     genelist=cellstr(genelist);
% end
lastwarn('')
save('input.mat','X','-v7.3');
writematrix(genelist,'g.txt');
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)    
    disp(warnId)
	if exist('./input.mat','file'), delete('./input.mat'); end
    disp('Writing data into input.txt...')
    sc_writefile('input.txt',X,genelist);
end
pkg.RunRcode('scriptnew.R');
if exist('output.h5','file')
    X=h5read('output.h5','/X');
elseif exist('output.txt','file')
    X=readmatrix('output.txt');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

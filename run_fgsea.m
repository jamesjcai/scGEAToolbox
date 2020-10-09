function [s]=run_fgsea(genelist,rmribo)
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end
if nargin<2, rmribo=false; end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_fgsea');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);


[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end



t=readtable('input_template.txt');


if exist('output.txt','file'), delete('output.txt'); end
genelist=upper(genelist);
a=-log(1-rand(length(genelist),1));
sortid=(1:length(genelist))';
% a=12+randn(length(genelist),1);
drdist=sort(a,'descend');
drdist=t.drdist(1:length(genelist));
T=table(sortid,genelist,drdist);
t=t(1:length(genelist),:);
t.genelist=T.genelist;
T=t;
if rmribo
    [gribo]=get_ribosomalgenes;
    i=~ismember(T.genelist,gribo);
    T=T(i,:);
end
writetable(T,'input.txt');
RunRcode('script.R');
pause(1);
if exist('output.txt','file')
    s=readtable('output.txt',"Delimiter",',');
else
    s=[];
end
% if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end
cd(oldpth);
end
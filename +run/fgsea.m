function [s]=fgsea(genelist,rmribo)

if nargin<2, rmribo=true; end
oldpth=pwd();
[isok,msg]=commoncheck_R('R_fgsea');
if ~isok, error(msg); end


if exist('output.txt','file'), delete('output.txt'); end


t=readtable('input_template.txt');
N=size(t,1);
if length(genelist)>N
    genelist=genelist(1:N);
end

genelist=upper(genelist);
%a=-log(1-rand(length(genelist),1));
sortid=(1:length(genelist))';
% a=12+randn(length(genelist),1);
% drdist=sort(a,'descend');

drdist=t.drdist(1:length(genelist));
T=table(sortid,genelist,drdist);
t=t(1:length(genelist),:);
t.genelist=T.genelist;
T=t;
if rmribo
    [gribo]=pkg.i_get_ribosomalgenes;
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
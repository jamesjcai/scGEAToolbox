function [s]=gseapy(genelist,drdistin)

if nargin<2, drdistin=[]; end
oldpth=pwd();
pw=fileparts(mfilename('fullpath'));
wkfold=fullfile(pw,'thirdparty','gseapy');
cd(wkfold);

if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end

t=readtable('input_template.txt');
N=size(t,1);
if length(genelist)>N
    genelist=genelist(1:N);
end

genelist=upper(genelist);
if ~isempty(drdistin) && length(genelist)==length(drdistin)
    drdist=drdistin;
    T=table(genelist,drdist);
else
    drdist=t.drdist(1:length(genelist));
    T=table(genelist,drdist);
end
writetable(T,'input.txt','WriteVariableNames',false);
% https://gseapy.readthedocs.io/en/latest/gseapy_example.html

%if exist('input.txt','file'), delete('input.txt'); end
%if exist('output.txt','file'), delete('output.txt'); end
%cd(oldpth);
end
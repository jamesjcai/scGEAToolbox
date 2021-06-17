function [A]=scode(X,t)

% if nargin<2, plotit=false; end
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end

oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_SCODE');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end



if size(t,2)==1
    t=[(1:length(t))' t];
end

%if ~exist('input.csv','file')
writematrix(X,'input1.txt');
writematrix(t,'input2.txt');
%end
RunRcode('script.R');
if exist('output.csv','file')
    A=readmatrix('output.csv');
    %G=digraph(A);    
else
    A=[];
    %G=[];
end

cd(oldpth);

% if plotit
%     if ~isempty(A)
%         plot(G);
%     end
% end
end
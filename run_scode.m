function [A]=run_scode(X,t)

% if nargin<2, plotit=false; end
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_SCODE');
cd(pth);

if size(t,2)==1
    t=[(1:length(t))' t];
end

%if ~exist('input.csv','file')
writematrix(X,'input1.txt');
writematrix(t,'input2.txt');
%end
RunRcode('script.R');
if exist('output.csv','file')
    A=csvread('output.csv',1,1);
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

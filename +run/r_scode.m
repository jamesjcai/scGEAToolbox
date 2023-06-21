function [A]=r_scode(X,t)

% if nargin<2, plotit=false; end
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end

oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_SCODE');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);
Rpath=getpref('scgeatoolbox','rexecutablepath');
[~,cmdout]=pkg.RunRcode('require.R',Rpath);
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end



if size(t,2)==1
    t=[(1:length(t))' t];
end

%if ~exist('input.csv','file')
if issparse(X), X=full(X); end
writematrix(X,'input1.txt');
writematrix(t,'input2.txt');
%end

pkg.RunRcode('script.R',Rpath);
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

function [A]=r_pcnet(X)

% if nargin<2, plotit=false; end
if isempty(pkg.FindRpath)
   error('Rscript.exe is not found. Use native matlab function SC_PCNET.m instead.');
end

oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_dna_differential_network_analysis');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

if exist('output.csv','file')
    delete('output.csv');
end

%if ~exist('input.csv','file')
if issparse(X), X=full(X); end
    writematrix(X,'input.csv');
%end

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if exist('output.csv','file')
    A=readmatrix('output.csv');
    %G=digraph(A);    
else
    A=[];
    %G=[];
end
if exist('input.csv','file')
    delete('input.csv');
end
cd(oldpth);

% if plotit
%     if ~isempty(A)
%         plot(G);
%     end
% end
end

function [A]=r_pcnet(X)

% if nargin<2, plotit=false; end
if isempty(FindRpath)
   error('Rscript.ext is not found. Use native matlab function SC_PCNET.m instead.');
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
    writematrix(X,'input.csv');
%end
RunRcode('script.R');
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
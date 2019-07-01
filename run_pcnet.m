function [A]=run_pcnet(X)

% if nargin<2, plotit=false; end
if isempty(FindRpath)
   error('Rscript.ext is not found. Use native matlab function SC_PCNET.m instead.');
end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_dna_differential_network_analysis');
cd(pth);

if exist('output.csv','file')
    delete('output.csv');
end

%if ~exist('input.csv','file')
    csvwrite('input.csv',X);
%end
RunRcode('script.R');
if exist('output.csv','file')
    A=csvread('output.csv',1,1);
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

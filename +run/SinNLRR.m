function [C]=SinNLRR(X,k)
% SinNLRR - 
%
% USAGE:
% >> % [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% load('example_data/example10xdata.mat');
% [C,s]=run_simlr(X,[],true);
% figure;
% scatter(s(:,1),s(:,2),20,C,'filled')

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','SinNLRR');
if ~(ismcc || isdeployed)
    addpath(pth);
end
if nargin<2 || isempty(k)
    k=fun_num_cluster(X);
    fprintf('k=%d\n',k);
end

% if nargin<3
%     donorm=false;
% end

% if donorm
%     % [X]=sc_norm(X);
%     % X=log10(X+1);
%     X=X./vecnorm(X);
% end
[C] = SinNLRRori(X,k);


end
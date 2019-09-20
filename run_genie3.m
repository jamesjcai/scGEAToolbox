function [A,s,G]=run_genie3(X,genelist,donorm,plotit)
%RUN_GENIE3
%
% USAGE:
% >>[X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >>[X,genelist]=sc_selectg(X,genelist,5,3);
% >>figure; [A,s,G]=run_genie3(X(1:8,:),genelist(1:8),false,true);

if nargin<2 || isempty(genelist) 
    genelist=string(1:size(X,1));
end
if nargin<3, donorm=false; end
if nargin<4, plotit=false; end

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/GENIE3');
addpath(pth);
pth=fullfile(pw1,'thirdparty/GENIE3/RT');
addpath(pth);

if donorm   
   X=sc_norm(X,'type','libsize');
end
% GENIE3 takes data (rows: cell, columns: genes)
data=X';

A=GENIE3(data);
if nargout==2
    s=get_link_list(A);
end
%%
if plotit || nargout==3
    s=get_link_list(A);
    G=digraph(s(:,1),s(:,2),s(:,3),genelist);
end
if plotit
    LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    LWidths(LWidths==0)=1e-5;
    % p=plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)
    % LWidths=ones(size(LWidths));
    try
        p=plot(G,'LineWidth',LWidths);
        % axis square
        p.MarkerSize = 7;
        p.Marker = 's';
        p.NodeColor = 'r';
    catch ME
        warning('Error in plotting.');
    end
end

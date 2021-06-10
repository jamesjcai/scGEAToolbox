function [Tp,Tn]=GSEA(genelist,generank,showplot,perm_nb,dbfile,sort_type)

if nargin<2, generank=1:length(genelist); end
if nargin<3, showplot=false; end
if nargin<4, perm_nb=1000; end
if nargin<5, dbfile='msigdb_v62_c5.mat'; end
if nargin<6, sort_type='ascend'; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','GSEA');
addpath(pth);

% set options
opts = default_GSEA_opts();
opts.show = showplot;      % if plot results
opts.save = false;         % if save results
opts.perm_nb = perm_nb;    % number of permutations
opts.sort_type=sort_type;
opts.GS_sort_type='NES_qval';

load(dbfile,'GeneSet','GeneSetName');
%  GeneSet=GeneSet(1:50);
%  GeneSetName=GeneSetName(1:50);

[res_pos,res_neg,res_descr] = MrGSEAPreranked(generank,genelist,...
                                     GeneSet,GeneSetName,'GSEA_result',opts);

% for k=1:size(res_pos,1)
%     res_pos{k,1}=GeneSetName{res_pos{k,1}};
% end
% for k=1:size(res_neg,1)
%     res_neg{k,1}=GeneSetName{res_neg{k,1}};
% end

% load GSEA_result
Tp=cell2table(res_pos);
Tp.Properties.VariableNames=res_descr;
Tp=sortrows(Tp,'NES_qval','ascend');
Tp.GSETNAME=GeneSetName(Tp.GSID);
Tp=movevars(Tp,'GSETNAME','before','GSID');

Tn=cell2table(res_neg);
Tn.Properties.VariableNames=res_descr;
Tn = sortrows(Tn,'NES_qval','ascend');
Tn.GSETNAME=GeneSetName(Tn.GSID);
Tn=movevars(Tn,'GSETNAME','before','GSID');


end
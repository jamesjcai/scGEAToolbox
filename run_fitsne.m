function s=run_fitsne(X,ndim,plotit,donorm,dolog1p,bygene,genelist)
%RUN_PHATE 
%
% PHATE is a data reduction method specifically designed for visualizing 
% **high** dimensional data in **low** dimensional spaces. Details on 
% this package can be found here https://github.com/KrishnaswamyLab/PHATE. 
% For help, visit https://krishnaswamylab.org/get-help. For a more in 
% depth discussion of the mathematics underlying PHATE, see the bioRxiv 
% paper here: https://www.biorxiv.org/content/early/2017/12/01/120378.
%
% USAGE:
% >>[X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >>[X,genelist]=sc_selectg(X,genelist);
% >>figure; s1=run_phate(X,3,true,true,genelist);  % view genes
% >>figure; s2=run_phate(X,3,false,true);          % view cells
% >>figure; scatter3(s2(:,1),s2(:,2),s2(:,3),10,'filled'); % using output S

if nargin<2, ndim=2; end
if nargin<3, plotit=false; end
if nargin<4, donorm=true; end
if nargin<5, dolog1p=true; end
if nargin<6, bygene=false; end   % when BYGENE=true, the matrix X will be transposed and the output will be tSNE for genes rather than cells.
if nargin<7, genelist=[]; end


pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/FIt-SNE');
addpath(pth);

if donorm
if bygene   
   data=sc_norm(X','type','libsize');
else
   data=sc_norm(X,'type','libsize');
end
end

% % The following transpose is necessary to make the input dim right.
data=data';

% sqrt transform
% data = sqrt(data);
if dolog1p
    data = log(data+1);
end

opts.no_dims=ndim;
opts.max_iter = 400;
s = fast_tsne(data, opts);

%%
if plotit
    if bygene
        [lgu, dropr, lgcv] = sc_stat(X, [], false);
        colorby='mean';
        switch colorby
            case 'mean'
                C=lgu;
            case 'dropout'
                C=dropr;
            case 'cv'
                C=lgcv;
        end
    else
        C=sum(X,1);
    end
    switch ndim
        case 2
            scatter(s(:,1), s(:,2), 10, C, 'filled');
            xlabel 'tSNE1'
            ylabel 'tSNE2'
            title 'tSNE'
        case 3            
            scatter3(s(:,1), s(:,2), s(:,3), 10, C, 'filled');
            xlabel 'tSNE1'
            ylabel 'tSNE2'
            zlabel 'tSNE3'
            title 'tSNE'

    end
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist};
    end    
end
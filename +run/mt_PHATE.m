function s=mt_PHATE(X,ndim,plotit,bygene,genelist)
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
% >>figure; s2=run_phate(X,3,true);          % view cells
% >>figure; scatter3(s2(:,1),s2(:,2),s2(:,3),10,'filled'); % using output S

if nargin<2, ndim=3; end
if nargin<3, plotit=false; end
if nargin<4, bygene=false; end
if nargin<5, genelist=[]; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'external','mt_PHATE');
if ~(ismcc || isdeployed)
    addpath(pth);
end
% gene_names=cellstr(gl123);  
% PHATE on data (rows: samples, columns: features)
% data=X';
% % library size normalization
% libsize = sum(data,2);
% data = bsxfun(@rdivide, data, libsize) * median(libsize);

if bygene   
   data=sc_norm(X','type','libsize');
else
   data=sc_norm(X,'type','libsize');
end

% The following transpose is necessary to make the input dim right.
data=data';

% sqrt transform
%data = log(data+1);
data = sqrt(data);


s = phate(data, 't', 20, 'ndim', ndim, 'k', 10);
%s = phate(data, 'ndim', ndim);
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
        C=sum(X,1);  % by library size
    end
    switch ndim
        case 2
            scatter(s(:,1), s(:,2), 10, C, 'filled');
            % colormap(jet)
            % set(gca,'xticklabel',[]);
            % set(gca,'yticklabel',[]);
            % axis tight
            xlabel 'PHATE1'
            ylabel 'PHATE2'
            title 'PHATE'
%             h = colorbar;
%             set(h,'xtick',1:5);
%             ylabel(h, 'time');
        case 3
            scatter3(s(:,1), s(:,2), s(:,3), 10, C, 'filled');
            % colormap(jet)
            % set(gca,'xticklabel',[]);
            % set(gca,'yticklabel',[]);
            % set(gca,'zticklabel',[]);
            % axis tight
            xlabel 'PHATE1'
            ylabel 'PHATE2'
            zlabel 'PHATE3'
            title 'PHATE 3D'
%           h = colorbar;
%           set(h,'xtick',1:5);
%           ylabel(h, 'time');
%           view([-170 15]);
    end
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist};
    end    
end
end
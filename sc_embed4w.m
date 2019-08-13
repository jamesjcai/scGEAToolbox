function [Y_PCA,Y_tSNE,Y_PHATE_2D,Y_PHATE_3D]=sc_embed4w(X,bygene,plotit,genelist,colorby)
if nargin<2, bygene=true; end
if nargin<3, plotit=true; end
if nargin<4, genelist=[]; end
if nargin<5, colorby='mean'; end

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/PHATE');
addpath(pth);

if bygene
   data=sc_norm(X','type','libsize');
else
   data=sc_norm(X,'type','libsize');
end
% sqrt transform
data = sqrt(data');
% data = log(data+0.1);

%% 
% PCA
tic;
Y_PCA = svdpca(data, 2, 'random');
toc

% size(Y_PCA)
% size(X)
% C=sum(X);
% size(C)
% scatter(Y_PCA(:,1), Y_PCA(:,2), 5, C, 'filled');
% return;


%     colormap(jet)
%     set(gca,'xticklabel',[]);
%     set(gca,'yticklabel',[]);
%     axis tight
%     xlabel 'PCA1'
%     ylabel 'PCA2'
%     title 'PCA'
% 
%     if ~isempty(genelist)
%         dt = datacursormode;
%         dt.UpdateFcn = {@i_myupdatefcn3,genelist,X};
%     end
%     return;

% PHATE 2D
% tic;
% Y_PHATE_2D = phate(data, 't', 20);
% toc;

% PHATE 3D
tic;
Y_PHATE_3D = phate(data, 't', 20, 'ndim', 3);
toc;

Y_PHATE_2D=Y_PHATE_3D(:,1:2);

% tSNE -- slow!!!
tic;
Y_tSNE = tsne(data,'Theta',0.5,'NumPCAComponents',100,'Verbose',2, 'perplexity', 20);
toc

%% plot combined
if plotit
    if bygene
        [lgu, dropr, lgcv] = sc_stat(X, [], false);
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
    % cells=zeros(size(Y_PCA,1),1);
    subplot(2,2,1);
    scatter(Y_PCA(:,1), Y_PCA(:,2), 5, C, 'filled');
%     colormap(jet)
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    axis tight
    xlabel 'PCA1'
    ylabel 'PCA2'
    title 'PCA'
%     h = colorbar;
%     set(h,'xtick',1:5);
%     ylabel(h, 'time');
%     box on
%     if ~isempty(genelist)
%         dt = datacursormode;
%         dt.UpdateFcn = {@i_myupdatefcn1,genelist};
%     end

    subplot(2,2,2);
    scatter(Y_tSNE(:,1), Y_tSNE(:,2), 5, C, 'filled');
%     colormap(jet)
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    axis tight
    xlabel 'tSNE1'
    ylabel 'tSNE2'
    title 'tSNE'
%     h = colorbar;
%     set(h,'xtick',1:5);
%     ylabel(h, 'time');
%     box on
%     if ~isempty(genelist)
%         dt = datacursormode;
%         dt.UpdateFcn = {@i_myupdatefcn1,genelist};
%     end

    subplot(2,2,3);
    scatter(Y_PHATE_2D(:,1), Y_PHATE_2D(:,2), 5, C, 'filled');
%     colormap(jet)
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    axis tight
    xlabel 'PHATE1'
    ylabel 'PHATE2'
    title 'PHATE 2D'
%     h = colorbar;
%     set(h,'xtick',1:5);
%     ylabel(h, 'time');
%     box on
%     if ~isempty(genelist)
%         dt = datacursormode;
%         dt.UpdateFcn = {@i_myupdatefcn1,genelist};
%     end

    subplot(2,2,4);
    scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 5, C, 'filled');
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    set(gca,'zticklabel',[]);
    axis tight
    xlabel 'PHATE1'
    ylabel 'PHATE2'
    zlabel 'PHATE3'
    title 'PHATE 3D'
%    view([-170 15]);
%     h = colorbar;
%     set(h,'xtick',1:5);
%     ylabel(h, 'time');
%     box on
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn3,genelist,X};
    end
end

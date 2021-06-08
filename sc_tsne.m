function s=sc_tsne(X,ndim,plotit,donorm,dolog1p,bygene,genelist)
%tSNE embedding of cells

%see also: RUN_PHATE, RUN_UMAP
% s_phate=run.PHATE(X,3,true);
% s_umap=run.UMAP(X,3);

if nargin<2, ndim=3; end
if nargin<3, plotit=false; end
if nargin<4, donorm=true; end
if nargin<5, dolog1p=false; end
if nargin<6, bygene=false; end   % when BYGENE=true, the matrix X will be transposed and the output will be tSNE for genes rather than cells.
if nargin<7, genelist=[]; end


pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'+run','thirdparty','PHATE');
addpath(pth);

if donorm
    if bygene   
       data=sc_norm(X','type','libsize');
    else
       data=sc_norm(X,'type','libsize');
    end
else
    if bygene
       data=X';
    else
       data=X;
    end
end

% The following transpose is necessary to make the input dim right.
data=data';

% sqrt transform
% data = sqrt(data);
if dolog1p
    data = log(data+1);
end
if issparse(data)
   data=full(data);
end

ncells=size(data,1);

if ncells>500
	data = svdpca(data, 50, 'random');
end

if ncells>5000
    s=tsne(data,'NumDimensions',ndim,'Algorithm','barneshut','NumPCAComponents',100);
else
    s=tsne(data,'NumDimensions',ndim);
end

%%
if plotit
    if bygene
        [lgu, dropr, lgcv] = sc_genestat(X, [], false);
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
        otherwise            
            scatter3(s(:,1), s(:,2), s(:,3), 10, C, 'filled');
            xlabel 'tSNE1'
            ylabel 'tSNE2'
            zlabel 'tSNE3'
            title 'tSNE 3D'
    end
    if bygene && ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist};
    end    
end
end

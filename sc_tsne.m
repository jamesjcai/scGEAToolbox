function s=sc_tsne(X,ndim,donorm,dolog1p)
%tSNE embedding of cells

%see also: RUN_PHATE, RUN_UMAP
% s_phate=run.PHATE(X,3,true);
% s_umap=run.UMAP(X,3);

if nargin<2, ndim=3; end
if nargin<3, donorm=true; end
if nargin<4, dolog1p=true; end

% if nargin<5, bygene=false; end   % when BYGENE=true, the matrix X will be transposed and the output will be tSNE for genes rather than cells.
% if nargin<6, genelist=[]; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'+run','thirdparty','PHATE');
if ~(ismcc || isdeployed), addpath(pth); end

%if bygene, X=X.'; end
if donorm
    X=sc_norm(X);
    disp('Library-size normalization...done.')
end

X(isnan(X))=0;     % empty cells will be retained to keep spots

if dolog1p
    X=log(X+1);
    disp('Log(x+1) transformation...done.')    
end

if issparse(X), X=full(X); end
[ngenes]=size(X,1);

% The following transpose is necessary to make the input dim right.
data=X.';
%if ncells>500
%	data = svdpca(data, 50, 'random');
%end
if ngenes>500
    s=tsne(data,'NumDimensions',ndim,...
        'Algorithm','barneshut','NumPCAComponents',50,...
        'Standardize',true);
else
    s=tsne(data,'NumDimensions',ndim,...
        'Algorithm','exact','NumPCAComponents',0,...
        'Standardize',true);
end

%%
% if plotit
%     if bygene
%         [lgu, dropr, lgcv] = sc_genestat(X, [], false);
%         colorby='mean';
%         switch colorby
%             case 'mean'
%                 C=lgu;
%             case 'dropout'
%                 C=dropr;
%             case 'cv'
%                 C=lgcv;
%         end
%     else
%         C=sum(X,1);
%     end
%     switch ndim
%         case 2
%             scatter(s(:,1), s(:,2), 10, C, 'filled');
%             xlabel 'tSNE1'
%             ylabel 'tSNE2'
%             title 'tSNE'
%         otherwise            
%             scatter3(s(:,1), s(:,2), s(:,3), 10, C, 'filled');
%             xlabel 'tSNE1'
%             ylabel 'tSNE2'
%             zlabel 'tSNE3'
%             title 'tSNE 3D'
%     end
%     if bygene && ~isempty(genelist)
%         dt = datacursormode;
%         dt.UpdateFcn = {@i_myupdatefcn1,genelist};
%     end
end

end

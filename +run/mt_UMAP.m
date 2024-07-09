function [s] = mt_UMAP(X, ndim)

%   addpath /Users/Stephen/umap
%   addpath /Users/Stephen/util
%   javaaddpath('/Users/Stephen/umap/umap.jar');


%if nargin<4, verbose=false; end
%if nargin<3, plotit=false; end
if nargin < 2, ndim = 3; end

pw1 = fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    pth1 = fullfile(pw1, 'external', 'mt_UMAP');
    addpath(pth1);
    %pth2 = fullfile(pw1, 'external', 'mt_UMAP', 'util');
    %addpath(pth2);

    %pth3 = fullfile(pw1, 'external', 'mt_UMAP44', 'umap.jar');
    %javaaddpath(pth3);
end

data = transpose(X);

% data=transpose(sc_transform(X));
% data=transpose(sc_norm(X));

ncells = size(data, 1);
if ncells > 500
    if ~(ismcc || isdeployed)
        pth = fullfile(pw1, 'external', 'mt_PHATE');
        addpath(pth);
    end
    data = svdpca(data, 50, 'random');
end

%if ~ispc
%    [s] = run_umap_lite(data, 'n_components', ndim, 'verbose', false);
%else
    [s] = run_umap_lite_super(data, ndim);
%end
%{

if nargout>1 || plotit
    %if verbose
    [s,~,c]=run_umap_main(data, ...
        'n_components',ndim);
    % else
    %     [s,~,c]=run_umap_main(data, ...
    %         'n_components',ndim, ...
    %         'verbose','none');
    % end
else

    %if ~(ismcc || isdeployed)
    %    if verbose
    [s]=run_umap_lite(data,'n_components',ndim);
    %     else
    %         [s]=run_umap_lite(data,'n_components',ndim,'verbose','none');
    %     end
    % else
    %      umap=UMAP;
    %      umap.method='MEX';
    %      umap.n_components=ndim;
    %      umap.min_dist=0.3;
    %      umap.distance_func='euclidean';
    %      s = umap.fit_transform(data);
    % end
end

% if plotit && ~isempty(s)
%     gui.i_gscatter3(s,c);
%     xlabel('UMAP 1')
%     ylabel('UMAP 2')
% end
end

%}

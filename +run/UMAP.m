function [s,c]=UMAP(X,ndim,plotit,verbose)

%   addpath /Users/Stephen/umap
%   addpath /Users/Stephen/util
%   javaaddpath('/Users/Stephen/umap/umap.jar');


if nargin<4, verbose=false; end
if nargin<3, plotit=false; end
if nargin<2, ndim=3; end

pw1=fileparts(mfilename('fullpath'));
pth1=fullfile(pw1,'thirdparty','umapFileExchange','umap');
pth2=fullfile(pw1,'thirdparty','umapFileExchange','util');
pth3=fullfile(pw1,'thirdparty','umapFileExchange','umap','umap.jar');
addpath(pth1);
addpath(pth2);
javaaddpath(pth3);


data=transpose(sc_transform(X));
% data=transpose(sc_norm(X));

ncells=size(data,1);
if ncells>500
    pth=fullfile(pw1,'thirdparty','PHATE');
    addpath(pth);
	data = svdpca(data, 50, 'random');
end

if nargout>1 || plotit
    if verbose
        [s,~,c]=run_umap_main(data,'n_components',ndim);
    else
        [s,~,c]=run_umap_main(data,'n_components',ndim,'verbose','none');
    end
else
    if verbose
        [s]=run_umap_main(data,'n_components',ndim);
    else
        [s]=run_umap_main(data,'n_components',ndim,'verbose','none');
    end    
end

if plotit && ~isempty(s)
    gui.i_gscatter3(s,c);
    xlabel('UMAP 1')
    ylabel('UMAP 2')
end
end


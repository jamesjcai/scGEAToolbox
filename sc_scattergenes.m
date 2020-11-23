function sc_scattergenes(X,genelist,plottype,dofit)
%Scatter plots for genes

if nargin<2, genelist=[]; end
if nargin<3, plottype='mean_cv'; end
if nargin<4, dofit=false; end
switch plottype
    case 'mean_cv'
        scatter_mean_vs_cv(X,genelist,dofit);
    case 'mean_dropr'
        scatter_mean_vs_dropr(X,genelist,3);
        % scatter_mean_vs_dropr(X,genelist,1);
    case 'meanlg_varlg'
        scatter_meanlg_vs_varlg(X,genelist,dofit);
end
function [Tct]=sc_celltypebrushed(X,genelist,s,brushedData,species,organ)

% USAGE:

% s=sc_tsne(X,3);
% figure; sc_cellscatter(s)
% % get brushedData
% [Tct]=sc_celltypesbrushed(X,genelist,s,brushedData)
if nargin<6, organ='all'; end
if nargin<5, species='human'; end

if islogical(brushedData)
    i=brushedData;
else
    [~,i]=ismember(brushedData,s,'rows');
end
Xi=X(:,i);
[Xi,gi]=sc_selectg(Xi,genelist);
[Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);


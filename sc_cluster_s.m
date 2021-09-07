function [c_clustid]=sc_cluster_s(s,k,varargin)
%sc_cluster_s - cluster cells using cell embeding s
%
%see also: sc_cluster_x

%if min(size(s))>3, error('S is coordinates of dimensional reduction.'); end

if nargin<2, k=6; end
p = inputParser;
defaultType = 'kmeans';
validTypes = {'kmeans','kmedoids','dbscan',...
    'spectclust','snndpc','mbkmeans'};
% 
checkType = @(x) any(validatestring(x,validTypes));

checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);

addRequired(p,'s',@isnumeric);
addRequired(p,'k',checkK);
addOptional(p,'type',defaultType,checkType)
addOptional(p,'plotit',false,@islogical)
parse(p,s,k,varargin{:})
plotit=p.Results.plotit;

switch p.Results.type
    case {'spectralcluster','spectclust'}
        c_clustid=spectralcluster(s,k);
    case 'kmeans'
        c_clustid=kmeans(s,k);
    case 'kmedoids'
        c_clustid=kmedoids(s,k);
    case 'dbscan'
        warning('In development. Needs parameters');
    case 'snndpc'
        c_clustid=sc_snndpc(s,k);
    case 'mbkmeans'
        [~,~,c_clustid]=pkg.mbkmeans(s,k);
end

if plotit
    gui.i_gscatter3(s,c_clustid);
    hold on
    for i=1:k
        si=s(c_clustid==i,:);
        si=mean(si);
        if size(s,2)==3
            text(si(:,1),si(:,2),si(:,3),sprintf('%d',i),...
                'fontsize',20,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
        else
            text(si(:,1),si(:,2),sprintf('%d',i),...
                'fontsize',20,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
        end
    end
end

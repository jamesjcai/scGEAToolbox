function [id]=sc_clustshow(s,k,varargin)

if nargin<2, k=3; end
p = inputParser;
defaultType = 'spectralcluster';
validTypes = {'spectralcluster','kmeans','kmedoids','dbscan'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'s',@isnumeric);
addRequired(p,'k',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,s,k,varargin{:})

switch p.Results.type
    case 'spectralcluster'
        id=spectralcluster(s,k);
    case 'kmeans'
        id=kmeans(s,k);
    case 'kmedoids'
        id=kmedoids(s,k);
    case 'dbscan'
        warning('In development. Needs parameters');        
end

if size(s,2)==3
    scatter3(s(:,1),s(:,2),s(:,3),10,id);
elseif size(s,2)==2
    scatter(s(:,1),s(:,2),10,id);
end
hold on
for i=1:k    
    si=s(id==i,:);
    si=mean(si);
    if size(s,2)==3
        text(si(:,1),si(:,2),si(:,3),sprintf('%d',i),'fontsize',20);
    else
        text(si(:,1),si(:,2),sprintf('%d',i),'fontsize',20);
    end
end


function [Tct]=local_celltypebrushed(X,genelist,s,...
                brushedData,species,organ,database)

if nargin<7, organ='panglaodb'; end
if nargin<6, organ='all'; end
if nargin<5, species='mouse'; end

if islogical(brushedData)
    i=brushedData;
else
    [~,i]=ismember(brushedData,s,'rows');
end
Xi=X(:,i);
[Xi,gi]=sc_selectg(Xi,genelist);
if strcmpi(database,'clustermole')
    %disp('Using clustermole marker database')
    [Tct]=sc_celltypecaller_new(Xi,gi,[],'species',species);
elseif strcmpi(database,'panglaodb')
    %disp('Using panglaodb marker database')
    [Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
end
end

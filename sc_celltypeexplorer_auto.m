function sc_celltypeexplorer_auto(X,genelist,s,varargin)

   p = inputParser;
   addOptional(p,'species',"mouse",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["human","mouse"]));
   parse(p,varargin{:});
   species=p.Results.species;
   
c=sc_clustshow(s);
k=max(c);    
figure;    
if size(s,2)==3
    scatter3(s(:,1),s(:,2),s(:,3),10,c);
elseif size(s,2)==2
    scatter(s(:,1),s(:,2),10,c);
end
hold on
for i=1:k
    ptsSelected=s(c==k,:);
    [Tct]=sc_celltypebrushed(X,genelist,s,ptsSelected,species);    
    si=s(c==i,:);
    si=mean(si);
    if size(s,2)==3
        text(si(:,1),si(:,2),si(:,3),sprintf('%s',Tct.C1_Cell_Type{1}),...
            'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
    else
        text(si(:,1),si(:,2),sprintf('%s',Tct.C1_Cell_Type{1}),...
            'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
    end
end
function sc_celltypeexplorer_auto(X,genelist,s,varargin)

   p = inputParser;
   addOptional(p,'species',"mouse",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["human","mouse"]));
   addOptional(p,'organ',"all",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["all","heart","immunesystem","brain","pancreas"]));
   addOptional(p,'method',"alona",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["alona","singler"]));
   parse(p,varargin{:});
   species=p.Results.species;
   organ=p.Results.organ;
   method=p.Results.method;
   
    c=sc_clustshow(s,6,'plotit',false);   
    %rng(1234)
    %c=sc_clustshow(s,6,'type','kmedoids','plotit',false);
    k=max(c);    
    % figure;    
    if size(s,2)==3
        scatter3(s(:,1),s(:,2),s(:,3),10,c);
    elseif size(s,2)==2
        scatter(s(:,1),s(:,2),10,c);
    end
    hold on
    %[T]=sc_celltypecaller(X,genelist,c);

for i=1:k
    %ptsSelected=s(c==i,:);
    %[Tct]=sc_celltypebrushed(X,genelist,s,ptsSelected,species);    
    Xi=X(:,c==i);
    [Xi,gi]=sc_selectg(Xi,genelist);    
    si=s(c==i,:);
    si=mean(si);
    
    switch lower(method)
        case 'alona'
            [Tct]=sc_celltypecaller(Xi,gi,[],'species',species,'organ',organ);
            ctxt=Tct.C1_Cell_Type{1};
        case 'singler'
            ctxt=run_singler(Xi,gi);
            ctxt=ctxt(mode(grp2idx(ctxt)));
    end    
    if size(s,2)==3
        text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
            'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
    else
        text(si(:,1),si(:,2),sprintf('%s',ctxt),...
            'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
    end
end
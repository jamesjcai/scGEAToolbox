function scatter_mean_vs_dropr(X,genelist,methodid)
if nargin<3, methodid=3; end
if nargin<2, genelist=[]; end

if issparse(X), X=full(X); end
dropr=1-sum(X>0,2)./size(X,2);
switch methodid
    case 1
        % seurte?
        X(X==0)=nan;
        ulg=nanmean(log(X),2);
    case 2
        % https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0805-z#MOESM1        
        % double exponential model fit.
        X(X==0)=nan;
        ulg=log(nanmean(X,2));
    case 3
        % M3Drop
        ulg=nanmean(X,2);
end

    if methodid~=3
        scatter(ulg,dropr);
        xlabel('Mean log(Expr)');
    else
        semilogx(ulg,dropr,'o');
        xlabel('Mean, log');
    end
    grid on
    ylabel('Dropout rate (% of zeros)');
    box on
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1,genelist};
    
    

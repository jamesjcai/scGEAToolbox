function [X,genelist]=sc_filterg(X,genelist,cutoff,verbose)
if nargin<2, genelist=[]; end
if nargin<3, cutoff=1; end
if nargin<4, verbose=false; end

if cutoff<1.0   % e.g., 0.9
    if verbose
        fprintf('Discard genes with poor expression values (more than %d%% zeros in all cells).\n',...
            100*cutoff); 
    end
    r=sum(X~=0,2)./size(X,2);
    i=r>=cutoff;
else 
    [u]=sum(X,2);
    i=u>=cutoff;
    if verbose
        fprintf('Discard genes with poor expression values (with less than %d reads among all cells).\n',...
            cutoff);
    end    
end
% We discard cells with poor gene expression values (more than 90% zeros in all cells)
% As the default, we filter all the genes with less than 5 reads among 99% of the samples
X=X(i,:);
if ~isempty(genelist), genelist=genelist(i,:); end
if verbose, fprintf('%d genes removed.\n',sum(~i)); end
end
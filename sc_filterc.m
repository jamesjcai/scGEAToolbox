function [X,i]=sc_filterc(X,cutoff,verbose)
if nargin<2, cutoff=1; end
if nargin<3, verbose=false; end
[s]=sum(X,1);
i=s>=cutoff;
X=X(:,i);
if verbose, fprintf('%d samples (cells) removed.\n',sum(~i)); end




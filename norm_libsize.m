function [X]=norm_libsize(X,factorn)
X=double(X);
lbsz=nansum(X);
if nargin<2, factorn=nanmedian(lbsz); end
X=(X./lbsz)*factorn;

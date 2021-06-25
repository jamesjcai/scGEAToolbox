function [X]=norm_libsize(X,factorn)
%Library size normalization
X=double(X);
lbsz=nansum(X);
% if nargin<2, factorn=nanmedian(lbsz); end
if nargin<2, factorn=10000; end
X=(X./lbsz)*factorn;
end
function [X]=norm_libsize(X,factorn)
%Library size normalization

%aka: namely cell depth normalization to the mean cell depth, followed by log1p (log1pPF)
% https://www.biorxiv.org/content/10.1101/2022.05.06.490859v1.full

X=double(X);
lbsz=sum(X,'omitnan');
if nargin<2, factorn=mean(lbsz,'omitnan'); end
% if nargin<2, factorn=median(lbsz,'omitnan'); end
% if nargin<2, factorn=10000; end
X=(X./lbsz)*factorn;
end

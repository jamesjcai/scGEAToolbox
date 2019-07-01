function [X]=sc_selectc(X,cutoff)
if nargin<2, cutoff=1000; end
libsz=sum(X);

if cutoff>1.0
    i=libsz>=cutoff;
else
    i=libsz>=quantile(libsz,cutoff);
end

% function [X]=sc_selectc(X,lwprct,upprct)
% if nargin<3, upprct=1.0; end
% if nargin<2, lwprct=0.15; end
% if upprct~=1
%     i=libsz>=quantile(libsz,lwprct) & libsz<=quantile(libsz,upprct);
% else
%     i=libsz>=quantile(libsz,lwprct);
% end

X=X(:,i);


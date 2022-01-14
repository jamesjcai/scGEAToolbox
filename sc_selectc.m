function [X,keptidx]=sc_selectc(X,libszcutoff,gnnumcutoff)
if nargin<3, gnnumcutoff=500; end
if nargin<2, libszcutoff=1000; end
libsz=sum(X,1);
gnnum=sum(X>0,1);

if libszcutoff>1.0
    keptidx = (libsz>=libszcutoff) & (gnnum>=gnnumcutoff);
else
    keptidx = (libsz>=quantile(libsz,libszcutoff)) & (gnnum>=gnnumcutoff);
end

% function [X]=sc_selectc(X,lwprct,upprct)
% if nargin<3, upprct=1.0; end
% if nargin<2, lwprct=0.15; end
% if upprct~=1
%     i=libsz>=quantile(libsz,lwprct) & libsz<=quantile(libsz,upprct);
% else
%     i=libsz>=quantile(libsz,lwprct);
% end

X=X(:,keptidx);

end
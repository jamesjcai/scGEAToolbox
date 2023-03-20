function [Y]=e_embedbyd(D,ndim,methodid)

if nargin<2, ndim=2; end
if nargin<3, methodid=1; end

IniY = pkg.randmds(D, ndim);
opt = statset('display','iter');

switch methodid
    case 1     % tsne by distance matrix
        Y = pkg.e_tsnebyd(D,ndim,IniY);
    case 2     % metric MDS
        Y = mdscale(D,ndim,'options',opt,'start',IniY, ...
            'Criterion','metricstress');
    case 3     % non-metric MDS
         Y = mdscale(D,ndim,'options',opt,'start',IniY, ...
             'Criterion','stress');
end

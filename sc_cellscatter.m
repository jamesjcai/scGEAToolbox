function sc_cellscatter(s,grp)
if nargin<2, grp=ones(size(s(:,1),1)); end
scatter3(s(:,1),s(:,2),s(:,3));

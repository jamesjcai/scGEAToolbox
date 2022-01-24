%Filename:      plot_data.m
%Last edit:     Feb 11 2019
%Description:   Generates scatter plots and saves them in a file
%Inputs:
%               -P:
%
%               a N x 2 array of N 2 dimensional points
%
%               -names:
%
%               a N x 1 cell with the names of the classes to plot
%
%               -filename:
%
%               name of the file where the plot will be saved
%
%               -indices:
%
%               a N x 1 array indicating the points that will be greyed out
%
%
%
% -------------------------------------------------------------------------

function plot_data(P, names, filename, indices)

hold off
aux=0;
if nargin>3
    scatter(P(not(indices),1), P(not(indices),2), '.', 'black');
    hold on
    gscatter(P(indices,1), P(indices,2), names(indices,1));
    aux=1;
else
    gscatter(P(:,1), P(:,2), names(:,1));
end
colormap('jet')
c=get(gca, 'Children');
for i=1:size(c,1)-aux
    c(i).MarkerSize=10;
end
[h,icons] = legend(gca);
h.Location='bestoutside';

axis off
set(h,'FontSize',16);

saveas(gcf, filename)

end

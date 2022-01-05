%Filename:      scGeneFit.m
%Last edit:     Feb 11 2019
%Description:   Implements scGeneFit to select markers from gene expression
%               data. The markers are chosen so that they preserve a given
%               hierarchical clustering structure
%Inputs:
%               -data: 
%     
%               a G x C array of gene expression data where G denotes 
%               the number of genes and C is the number of cells.
% 
%               -labels: 
% 
%               a 1 x L cell where L is the number of levels of the 
%               hierarchy. labels{i} is a 1 x C array consisting of the
%               the labels of each cell at that point of the hierarchy
%               labels{i} has entries ranging from 1 to Li where Li is
%               the number of clusters at that level of the hierarchy 
%               Example: L={[1,1,1,1,2,2,2,2], [1,1,2,2,1,1,2,2]}
%
%               -num_markers: 
%     
%               number of markers to select
%
%               -opt: 
%     
%               struct of optional parameters:
%               opt.samples_fraction (default 1): 
%                   fraction of data to get the constraints from.
%                   if smaller than 1 the data is selected at random.
%               opt.constraints_neighbors (default 4)
%                   number of nearest neighbors to select constraints from
%               opt.hinge_scale (default 0.1)
%                   hinge parameter for linear program (see paper)
% 
%Outputs:
%               -markers: 
% 
%               a G x 1 binary array. G(t)=true implies that t is a selected marker
%
%               -projected_data: 
% 
%               a num_markers x C array with data projected to the selected markers 
% 
%Documentation:
% 
% -------------------------------------------------------------------------
function [markers, projected_data]=scGeneFit_centers(data, labels, num_markers, opt)


%number of levels of the hierarchy
num_levels=size(labels,2);

%Unpack optional arguments and set default
samples_fraction= ones(1, num_levels);
constraints_neighbors=4;
hinge_scale=0.1;
if nargin>3
    if isfield(opt, 'samples_fraction')
        samples_fraction=opt.samples_fraction;
    end
    if isfield(opt, 'constraints_neighbors')
        constraints_neighbors=opt.constraints_neighbors;
    end
    if isfield(opt, 'hinge_scale')
        hinge_scale=opt.hinge_scale;
    end
end
tic
rng default %for reproducibility

%Select constraints from the top level 
[Delta,smallest,lambda]=select_center_constraints(data, labels{1}, samples_fraction(1), 95);

%Select constraints from other levels
for level=2:num_levels

for i=1:max(labels{level-1})
[indices,sublabels]=subcluster_by_label(labels{level-1}, i, labels{level});
[D, s, lam]=select_center_constraints(data(:,indices), sublabels, samples_fraction(level), 95);
Delta=horzcat(Delta,D);
lambda=vertcat(lambda,lam);
if s<smallest
    smallest=s;
end
end
end

toc
disp('Size of the linear program constraints: ')
size(Delta)
M=sqz_hinge(Delta, hinge_scale*smallest, num_markers, lambda);
[~,aux]=maxk(M,num_markers);
markers=M*0;
markers(aux)=1;
markers=markers>0;
projected_data=data(markers,:);
end
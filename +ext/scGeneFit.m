function [markerlist]=scGeneFit(X,g,c,num_markers)
%Genetic marker selection using linear programming
%
% https://github.com/solevillar/scGeneFit
% https://www.nature.com/articles/s41467-021-21453-4

if nargin<4
    num_markers=length(unique(c))*5;
end
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','scGeneFit');
addpath(pth);
opt.samples_fraction=0.15;
opt.constraints_neighbors=3;
opt.hinge_scale=7;

[markers,proj_data]=scGeneFit_centers(X,{c}, num_markers, opt);
[thisc,~]=grp2idx(c);
z=[];
for kk=1:max(thisc)
    z=[z mean(proj_data(:,thisc==kk),2)];
end
[~,idx]=max(z,[],2);
markerlist=cell(max(thisc),1);
gm=g(markers);

for kk=1:max(thisc)     
    markerlist{kk}=gm(idx==kk);
end
end

%{
    [c,idx]=sort(c);
    proj_data=proj_data(:,idx);
    % proj_data=proj_data(:,1:5:end);
    proj_data=log(sc_norm(proj_data)+1);
    % T=array2table(proj_data);
    % T.Properties.RowNames=g(markers);
    z=[];
    uc=unique(c);
    for kk=1:numel(uc)
        z=[z mean(proj_data(:,c==uc(kk)),2)];
    end
    [~,idx]=max(z,[],2);
    [~,idx]=sort(idx);
    proj_data=proj_data(idx,:);    
    gm=matlab.lang.makeUniqueStrings(g(markers));
    gm=gm(idx);
    proj_data=proj_data(:,1:5:end);
    
    h=heatmap(num2cell(1:size(proj_data,2)),...
        cellstr(gm),proj_data);
    
    h.Title = 'Marker Gene Heatmap';
    h.XLabel = 'Cell Type';
    h.YLabel = 'Marker Gene';
    h.Colormap = parula;
    h.GridVisible = 'off';
    % sorty(h);
%}
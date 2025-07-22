function [H, clrs, scheme]=scatterByColorScheme(...
    scheme, vector, ax, xyz, sz, mrkr)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if ~isempty(scheme)
    try
        J=feval(scheme);
    catch
        scheme=[];
    end
end
if isempty(scheme)
    [J, scheme]=SuhHeatMap.ColorScheme;
end
nJ=size(J,1);
vector(vector<0)=0;
vector(isnan(vector))=0;
good=sum(vector>=0 & vector <=1);
if good/length(vector)<.9
    vector=MatBasics.NormalizeByStd(vector);
end
clrIdxs=ceil(vector*nJ);
clrIdxs(clrIdxs==0)=1;
clrs=J(clrIdxs,:); 
if nargin>2
    if nargin<6
        mrkr='o';
        if nargin<5
            sz=22;
        end
    end
    if size(xyz, 2)==2
        H=scatter(ax, xyz(:,1), xyz(:,2), sz, clrs, mrkr, 'filled');
    else
        H=scatter3(ax, xyz(:,1), xyz(:, 2),xyz(:,3), sz, clrs, ...
            mrkr, 'filled');
    end
else
    H=[];
end
end
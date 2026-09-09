function [T, X1, X2, g, xyz1, xyz2, ...
        px1, py1, pz1, px2, py2, pz2] = e_dvanalysis_brennecke(sce1, sce2, cL1, cL2)
if nargin < 3, cL1 = {'1'}; end
if nargin < 4, cL2 = {'2'}; end
[T, X1, X2, g, xyz1, xyz2, px1, py1, pz1, px2, py2, pz2] = ...
    sc_dvg(sce1, sce2, cL1, cL2, 'brennecke');
end

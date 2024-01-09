function [sce] = sc_csubtypeanno(sce, cell_type_target, speciestag)
    if nargin < 3 || isempty(speciestag), speciestag='human'; end
    pw1 = fileparts(mfilename('fullpath'));
    pth = fullfile(pw1, 'resources', 'cellsubtypes.xlsx');
end
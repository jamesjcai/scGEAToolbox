function [primarymarkerstr] = e_primarymarkers(cell_type_target, speciestag)
%E_PRIMARYMARKERS Positive markers of one primary cell type, by species.
%   s = pkg.e_primarymarkers(cell_type_target) returns the comma-separated
%   marker string PanglaoDB lists for that primary cell type in human.
%
%   s = pkg.e_primarymarkers(cell_type_target, speciestag) selects
%   'human' (default) or 'mouse'.
%
%   The two marker sets are not the same list: over the primary types the
%   subtype table covers, the shared fraction runs from 0.73 to 1.00, and
%   4 primary types exist only in the human file against 8 only in the
%   mouse one. Picking the wrong file therefore changes which markers a
%   subtype is scored on rather than merely renaming them.
%
%   This was a local function inside SC_CSUBTYPEANNO that ignored the
%   species entirely and always loaded marker_hs.mat. It is out here so
%   the species choice can be tested without re-embedding a dataset.
%
%   See also SC_CSUBTYPEANNO, PKG.E_MARKERWEIGHT, PKG.E_DETERMINECELLTYPE.

if nargin < 2 || isempty(speciestag)
    speciestag = 'human';
end
speciestag = validatestring(lower(char(speciestag)), {'human', 'mouse'}, ...
    'pkg.e_primarymarkers', 'speciestag', 2);

switch speciestag
    case 'mouse'
        fname = 'marker_mm.mat';
    otherwise
        fname = 'marker_hs.mat';
end

pw1 = fileparts(fileparts(mfilename('fullpath')));   % .../+pkg -> root
pth1 = fullfile(pw1, 'external', 'fun_alona_panglaodb', fname);
load(pth1, 'Tm');

idx = upper(string(Tm.Var1)) == upper(string(cell_type_target));
if ~any(idx)
    % Naming the species matters: the same target can be valid for one
    % file and absent from the other.
    error('pkg:e_primarymarkers:UnknownPrimaryType', ...
        '"%s" is not a primary cell type in the %s marker set (%s).', ...
        string(cell_type_target), speciestag, fname);
end

primarymarkerstr = Tm.Var2{find(idx, 1)};
primarymarkerstr = strtrim(primarymarkerstr);
primarymarkerstr = erase(primarymarkerstr, " ");
primarymarkerstr = strip(primarymarkerstr, 'right', ',');
end

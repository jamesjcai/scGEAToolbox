function roles = i_h5adobsroles(names)
%I_H5ADOBSROLES  Which .h5ad obs columns the importer consumes, by index.
%
%   roles = pkg.i_h5adobsroles(names)
%
%   names  the obs column names, as PKG.I_READH5ADOBS returns them
%   roles  struct with fields "batch" and "celltype", each the index into
%          NAMES of the column playing that role, or 0 when absent
%
% Depositors name these columns differently, so each role has a list of
% aliases and names are compared ignoring case and any non-alphanumeric
% character: "CellType", "cell_type" and "cell type" all match the same
% entry, and only genuinely different words need listing.
%
% Both sides of the import ask here. SC_READH5ADFILE reads these columns
% out as its batch and cell type outputs, and PKG.I_ADDH5ADOBSATTRIBS
% leaves them out of the cell attributes it keeps, so a column cannot end
% up on the object twice under two spellings. Matching by name rather than
% by value is what makes that hold: a caller that replaces unlabeled cells
% with "undetermined" changes the values, so comparing them would no longer
% recognize the column it just consumed.
%
% See also SC_READH5ADFILE, PKG.I_READH5ADOBS, PKG.I_ADDH5ADOBSATTRIBS

arguments
    names string
end

aliases = struct( ...
    'batch',    {["BatchID", "batch", "sample_id", "donor_id", "orig.ident"]}, ...
    'celltype', {["CellType", "cell_types", "cell_ontology_class", ...
                  "cell_annotation", "annotation"]});

roles = struct('batch', 0, 'celltype', 0);
if isempty(names), return; end

simplify = @(s) lower(regexprep(s, '[^a-zA-Z0-9]', ''));
simplified = simplify(names(:)');

for role = string(fieldnames(aliases))'
    for candidate = aliases.(role)
        hit = find(simplified == simplify(candidate), 1);
        if ~isempty(hit)
            roles.(role) = hit;
            break;
        end
    end
end
end

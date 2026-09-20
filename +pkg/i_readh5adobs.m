function [names, values] = i_readh5adobs(filenm)
%I_READH5ADOBS  Read and decode every column of an .h5ad /obs table.
%
%   [names, values] = pkg.i_readh5adobs(filenm)
%
%   names   1xN string, the obs column names as the file spells them
%   values  1xN cell, each one value per cell
%
% The DataFrame index is left out: it holds the barcodes, which are read
% separately as the cell IDs, and it is not a column of the table.
%
% AnnData stores a column in one of three shapes, all handled here:
%   * a plain dataset, numeric or string
%   * a categorical group, "categories" plus 0-based "codes"
%   * a nullable group, "values" plus a logical "mask"
% A column in any other shape is skipped rather than guessed at.
%
% This is the one place the obs encodings are decoded. SC_READH5ADFILE
% picks its batch and cell-type columns out of what this returns, and
% PKG.I_ADDH5ADOBSATTRIBS keeps the rest as cell attributes, so a fix to
% the decoding reaches both.
%
% See also SC_READH5ADFILE, PKG.I_ADDH5ADOBSATTRIBS

arguments
    filenm (1,1) string
end

names = string.empty(1, 0);
values = {};

try
    info = h5info(filenm, '/obs');
catch
    return;   % no obs table; nothing to read
end

% The index is named by the group's "_index" attribute, defaulting to the
% conventional "_index" when the attribute is absent.
try
    indexName = string(h5readatt(filenm, '/obs', '_index'));
catch
    indexName = "_index";
end

% h5info gives leaf names for datasets but full paths for groups
members = string.empty(1, 0);
if ~isempty(info.Datasets)
    members = [members, string({info.Datasets.Name})];
end
if ~isempty(info.Groups)
    members = [members, regexprep(string({info.Groups.Name}), '^.*/', '')];
end

for k = 1:numel(members)
    name = members(k);
    if name == indexName, continue; end
    v = i_decodecolumn(filenm, "/obs/" + name);
    if isempty(v), continue; end
    names(end+1) = name; %#ok<AGROW>
    values{end+1} = v;   %#ok<AGROW>
end
end


function v = i_decodecolumn(h5file, path)
%I_DECODECOLUMN One obs column, whichever of the three shapes it is stored in.

% Shape 1: a plain dataset
try
    raw = h5read(h5file, char(path));
    if ischar(raw) || iscellstr(raw) || isstring(raw)
        raw = string(raw);
    end
    v = raw(:);
    return;
catch
    % Not a dataset; it is a group, so fall through
end

% Shape 2: categorical, categories plus 0-based codes
try
    codes = double(h5read(h5file, char(path + "/codes")));
    cats = string(h5read(h5file, char(path + "/categories")));
    cats = cats(:);

    % AnnData writes an unlabeled cell as an empty category name, which is
    % not a legal MATLAB category. Leaving those names out of the value set
    % makes the affected cells <undefined>, which is what they mean. A code
    % of -1 marks a missing value and falls outside the set the same way.
    named = strlength(strtrim(cats)) > 0;
    v = categorical(codes(:)+1, find(named), cats(named));
    return;
catch
    % Not categorical either
end

% Shape 3: nullable, values plus a logical mask marking the missing ones
try
    vals = h5read(h5file, char(path + "/values"));
    vals = vals(:);
    try
        mask = logical(h5read(h5file, char(path + "/mask")));
        if isnumeric(vals)
            vals = double(vals);
            vals(mask(:)) = NaN;
        elseif isstring(vals)
            vals(mask(:)) = missing;
        end
    catch
        % No mask; take the values as they stand
    end
    v = vals;
catch
    v = [];   % an encoding this does not know; the caller skips it
end
end

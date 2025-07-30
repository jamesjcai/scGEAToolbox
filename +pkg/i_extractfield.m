function a = i_extractfield(S, fieldname)
%EXTRACTFIELD_CUSTOM  Extract field values from struct array
%   a = EXTRACTFIELD_CUSTOM(S, fieldname) behaves like MATLAB's
%   extractfield, returning a 1×n numeric array if all field values
%   are numeric and uniform type, otherwise a 1×n cell array.
%
%   Inputs:
%     S         — structure array (any shape)
%     fieldname — string scalar or character vector naming the field
%
%   Output:
%     a — 1×numel(S) array or cell array of values

    validateattributes(S, {'struct'}, {'vector'}, mfilename, 'S', 1);
    validateattributes(fieldname, {'char', 'string'}, {'scalartext'}, mfilename, 'fieldname', 2);

    n = numel(S);
    out = cell(1, n);
    for k = 1:n
        out{k} = S(k).(fieldname);
    end

    % Check for uniform numeric type
    isAllNumeric = all(cellfun(@(x) isnumeric(x) && isscalar(x), out));
    if isAllNumeric
        try
            a = [out{:}];  % concatenate scalar numerics into vector
        catch
            a = out; % fallback to cell
        end
    else
        a = out;  % character vectors or mixed types → cell output
    end
end

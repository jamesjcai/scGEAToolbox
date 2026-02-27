function obj = setAnnotation(obj, fieldName, values)
    % Ensure valid name
    if ~isvarname(fieldName)
        error("Invalid annotation field name: %s", fieldName);
    end

    % Expand scalar
    if isscalar(values)
        values = repmat(string(values), obj.NumCells, 1);
    elseif numel(values) ~= obj.NumCells
        error("Annotation length mismatch: expected %d entries, got %d", ...
            obj.NumCells, numel(values));
    else
        values = string(values);
    end

    % Add/overwrite column in table
    obj.cellAnn.(fieldName) = values;
end
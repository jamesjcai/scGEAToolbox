function obj = subsetCellsByAnnotation(obj, fieldName, valuesToKeep)
    if ~ismember(fieldName, obj.cellAnn.Properties.VariableNames)
        error("Annotation field '%s' not found.", fieldName);
    end

    ann = obj.cellAnn.(fieldName);
    mask = ismember(ann, string(valuesToKeep));

    if ~any(mask)
        warning("No cells matched %s in field '%s'", strjoin(string(valuesToKeep), ","), fieldName);
    end

    obj = obj.subsetCells(mask);
end
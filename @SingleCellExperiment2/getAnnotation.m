function values = getAnnotation(obj, fieldName)
    if ismember(fieldName, obj.cellAnn.Properties.VariableNames)
        values = obj.cellAnn.(fieldName);
    else
        error("Annotation field '%s' not found.", fieldName);
    end
end
function uniqueNames = i_fastuniquenames(numFeatures)
% from metafeatures.m
    numcell = cellfun(@num2str,num2cell(1:numFeatures,1),'UniformOutput',false)';
    uniqueNames = strcat('Var',numcell);
end
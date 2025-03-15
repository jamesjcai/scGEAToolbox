function [T, Tnt] = in_DETableProcess(T, cL1, cL2)
    try
        T = sortrows(T, 'p_val_adj', 'ascend');
        T = sortrows(T, 'pct_1', 'ascend');
        T = sortrows(T, 'pct_2', 'descend');
        T = sortrows(T, 'avg_log2FC', 'ascend');
    catch ME
        warning(ME.message);
    end
    
    try
        if contains(T.Properties.VariableNames{5}, 'avg_1')
            T.Properties.VariableNames{5} = sprintf('%s_%s', ...
                T.Properties.VariableNames{5}, ...
                matlab.lang.makeValidName(string(cL1)));
        end
    
        if contains(T.Properties.VariableNames{6}, 'avg_2')
            T.Properties.VariableNames{6} = sprintf('%s_%s', ...
                T.Properties.VariableNames{6}, ...
                matlab.lang.makeValidName(string(cL2)));
        end
    
        if contains(T.Properties.VariableNames{7}, 'pct_1')
            T.Properties.VariableNames{7} = sprintf('%s_%s', ...
                T.Properties.VariableNames{7}, ...
                matlab.lang.makeValidName(string(cL1)));
        end
    
        if contains(T.Properties.VariableNames{8}, 'pct_2')
            T.Properties.VariableNames{8} = sprintf('%s_%s', ...
                T.Properties.VariableNames{8}, ...
                matlab.lang.makeValidName(string(cL2)));
        end
    catch ME
        warning(ME.message);
    end
        variables = T.Properties.VariableNames;
        for k = 1:length(variables)
            xx = T.(variables{k});
            if isnumeric(xx) && any(isinf(xx))
                xx(isinf(xx) & xx > 0) = 1e99;
                xx(isinf(xx) & xx < 0) = -1e99;
                T.(variables{k}) = xx;
            end
        end

        if nargout>1

            Item = T.Properties.VariableNames';
            Item = [Item; {'# of cells in sample 1';'# of cells in sample 2'}];
            Description = {'gene name';'p-value';...
                'log2 fold change between average expression';...
                'absolute value of log2 fold change';...
                'average expression in sample 1';...
                'average expression in sample 2';...
                'percentage of cells expressing the gene in sample 1';...
                'percentage of cells expressing the gene in sample 2';...
                'adjusted p-value'; 'test statistic'; ...
                 sprintf('%d',sum(i1&idx)); sprintf('%d',sum(i2&idx))};
            if length(Item) == length(Description)
                Tnt = table(Item, Description);
            else
                Tnt = table(Item);
            end
        end

end



%{

    function [T] = in_DETableProcess(T, cL1, cL2)
        try
            T = sortrows(T, 'p_val_adj', 'ascend');
            T = sortrows(T, 'pct_1', 'ascend');
            T = sortrows(T, 'pct_2', 'descend');
            T = sortrows(T, 'avg_log2FC', 'ascend');
        catch ME
            warning(ME.message);
        end
        
        try
            if contains(T.Properties.VariableNames{5}, 'avg_1')
                T.Properties.VariableNames{5} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{5}, ...
                    matlab.lang.makeValidName(string(cL1)));
            end
        
            if contains(T.Properties.VariableNames{6}, 'avg_2')
                T.Properties.VariableNames{6} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{6}, ...
                    matlab.lang.makeValidName(string(cL2)));
            end
        
            if contains(T.Properties.VariableNames{7}, 'pct_1')
                T.Properties.VariableNames{7} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{7}, ...
                    matlab.lang.makeValidName(string(cL1)));
            end
        
            if contains(T.Properties.VariableNames{8}, 'pct_2')
                T.Properties.VariableNames{8} = sprintf('%s_%s', ...
                    T.Properties.VariableNames{8}, ...
                    matlab.lang.makeValidName(string(cL2)));
            end
        catch ME
            warning(ME.message);
        end
            variables = T.Properties.VariableNames;
            for k = 1:length(variables)
                xx = T.(variables{k});
                if isnumeric(xx) && any(isinf(xx))
                    xx(isinf(xx) & xx > 0) = 1e99;
                    xx(isinf(xx) & xx < 0) = -1e99;
                    T.(variables{k}) = xx;
                end
            end
    end
%}
function [T, Tnt] = in_DVTableProcess(T, cL1, cL2)

        if nargout>1

            Item = T.Properties.VariableNames';
            Item = [Item; {'# of cells in sample 1';'# of cells in sample 2'}];
            
            Description = {'gene name';'log mean in sample 1';...
                'log CV in sample 1'; 'dropout rate in sample 1';...
                'distance to curve 1';'p-value of distance in sample 1';...
                'FDR of distance in sample 1';'log mean in sample 2';...
                'log CV in sample 2'; 'dropout rate in sample 2';...
                'distance to curve 2'; 'p-value of distance in sample 2';...
                'FDR of distance in sample 2'; 'Difference in distances';...
                'Sign of difference';'p-value of DV test';...
                sprintf('%d',sce1.NumCells); sprintf('%d',sce2.NumCells)};
            if length(Item) == length(Description)
                Tnt = table(Item, Description);
            else
                assignin("base","Item", Item);
                assignin("base","Description", Description);
                Tnt = table(Item);
                warning('Variables must have the same number of rows.');
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
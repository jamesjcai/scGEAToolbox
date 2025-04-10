function LassoAnalysisApp(X, y, varNames, parentfig)
    if nargin<4, parentfig = []; end
    % LassoAnalysisApp - LASSO Regression Analysis Tool with GUI
    % 
    % Syntax:
    %   LassoAnalysisApp(X, y)
    %   LassoAnalysisApp(X, y, varNames)
    % 
    % Inputs:
    %   X - Matrix of predictor variables (observations × variables)
    %   y - Vector of response variable
    %   varNames - (Optional) Cell array of variable names for predictors
    
    % Input validation
    if nargin < 2
        error('Both X (predictors) and y (response) must be provided');
    end
    
    % Validate dimensions
    [n_obs, n_vars] = size(X);
    if length(y) ~= n_obs
        error('Number of observations in X and y must match');
    end
    
    % Set variable names if not provided
    if nargin < 3 || isempty(varNames)
        varNames = cell(1, n_vars);
        for i = 1:n_vars
            varNames{i} = sprintf('Var%d', i);
        end
    elseif length(varNames) ~= n_vars
        error('Number of variable names must match number of predictor variables');
    end
    
    % Create the UI figure
    fig = uifigure('Name', 'LASSO Regression Analysis Tool', ...
        'Position', [0 0 600 600], 'Visible','off');

            try
               if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
                    [px_new] = gui.i_getchildpos(parentfig, fig);
                    if ~isempty(px_new)
                        movegui(fig, px_new);
                    else
                        movegui(fig, 'center');
                    end
                else
                    movegui(fig, 'center');
                end
            catch
                movegui(fig, 'center');
            end
    
    
    % Create the main layout - now with just 2 panels (parameters and results)
    mainLayout = uigridlayout(fig, [2 1]);
    mainLayout.RowHeight = {'fit', '1x'};
    
    % Create panels - removed data summary panel
    paramPanel = uipanel(mainLayout);
    paramPanel.Title = 'LASSO Parameters';
    paramPanel.Layout.Row = 1;
    
    resultsPanel = uipanel(mainLayout);
    resultsPanel.Title = 'Results';
    resultsPanel.Layout.Row = 2;
    
    % Parameter Panel Components
    paramGrid = uigridlayout(paramPanel, [9 2]); % Added row for cutoff control
    paramGrid.RowHeight = repmat({'fit'}, 1, 9);
    paramGrid.ColumnWidth = {'fit', '1x'};
    
    % Alpha parameter (LASSO = 1, Ridge = 0, Elastic Net = between)
    alphaTextLabel = uilabel(paramGrid);
    alphaTextLabel.Text = 'Alpha (0-1):';
    alphaTextLabel.Layout.Row = 1;
    alphaTextLabel.Layout.Column = 1;
    
    alphaDescLabel = uilabel(paramGrid);
    alphaDescLabel.Text = '1 = LASSO, 0 = Ridge';
    alphaDescLabel.Layout.Row = 1;
    alphaDescLabel.Layout.Column = 2;
    
    alphaSlider = uislider(paramGrid);
    alphaSlider.Limits = [0 1];
    alphaSlider.Value = 1;
    alphaSlider.Layout.Row = 2;
    alphaSlider.Layout.Column = [1 2];
    
    % Lambda selection method
    lambdaMethodLabel = uilabel(paramGrid);
    lambdaMethodLabel.Text = 'Lambda Selection:';
    lambdaMethodLabel.Layout.Row = 3;
    lambdaMethodLabel.Layout.Column = 1;
    
    lambdaMethodDropdown = uidropdown(paramGrid);
    lambdaMethodDropdown.Items = {'CV', 'MSE', 'Manual'};
    lambdaMethodDropdown.Value = 'CV';
    lambdaMethodDropdown.Layout.Row = 3;
    lambdaMethodDropdown.Layout.Column = 2;
    lambdaMethodDropdown.ValueChangedFcn = @updateLambdaControls;
    
    % Number of CV folds
    cvFoldsLabel = uilabel(paramGrid);
    cvFoldsLabel.Text = 'CV Folds:';
    cvFoldsLabel.Layout.Row = 4;
    cvFoldsLabel.Layout.Column = 1;
    
    cvFoldsSpinner = uispinner(paramGrid);
    cvFoldsSpinner.Value = 10;
    cvFoldsSpinner.Limits = [3 20];
    cvFoldsSpinner.Layout.Row = 4;
    cvFoldsSpinner.Layout.Column = 2;
    
    % Manual lambda value
    lambdaValueLabel = uilabel(paramGrid);
    lambdaValueLabel.Text = 'Lambda Value:';
    lambdaValueLabel.Layout.Row = 5;
    lambdaValueLabel.Layout.Column = 1;
    
    lambdaEdit = uieditfield(paramGrid, 'numeric');
    lambdaEdit.Value = 0.01;
    lambdaEdit.Layout.Row = 5;
    lambdaEdit.Layout.Column = 2;
    lambdaEdit.Enable = 'off';
    
    % Number of lambda values for path
    numLambdaLabel = uilabel(paramGrid);
    numLambdaLabel.Text = 'Num Lambda Values:';
    numLambdaLabel.Layout.Row = 6;
    numLambdaLabel.Layout.Column = 1;
    
    numLambdaSpinner = uispinner(paramGrid);
    numLambdaSpinner.Value = 100;
    numLambdaSpinner.Limits = [10 200];
    numLambdaSpinner.Layout.Row = 6;
    numLambdaSpinner.Layout.Column = 2;
    
    % Standardize data option
    standardizeLabel = uilabel(paramGrid);
    standardizeLabel.Text = 'Standardize:';
    standardizeLabel.Layout.Row = 7;
    standardizeLabel.Layout.Column = 1;
    
    standardizeCheck = uicheckbox(paramGrid);
    standardizeCheck.Value = true;
    standardizeCheck.Layout.Row = 7;
    standardizeCheck.Layout.Column = 2;
    
    % Variable cutoff selection - NEW
    cutoffLabel = uilabel(paramGrid);
    cutoffLabel.Text = 'Variable Selection:';
    cutoffLabel.Layout.Row = 8;
    cutoffLabel.Layout.Column = 1;
    
    cutoffPanel = uipanel(paramGrid, 'BorderType', 'none');
    cutoffPanel.Layout.Row = 8;
    cutoffPanel.Layout.Column = 2;
    
    cutoffGrid = uigridlayout(cutoffPanel, [2 2]);
    cutoffGrid.RowHeight = {'fit', 'fit'};
    cutoffGrid.ColumnWidth = {'fit', '1x'};
    cutoffGrid.Padding = [0 0 0 0];
    
    % Cutoff method
    cutoffMethodLabel = uilabel(cutoffGrid);
    cutoffMethodLabel.Text = 'Method:';
    cutoffMethodLabel.Layout.Row = 1;
    cutoffMethodLabel.Layout.Column = 1;
    
    cutoffMethodDropdown = uidropdown(cutoffGrid);
    cutoffMethodDropdown.Items = {'Coefficient Threshold', 'Top N Variables'};
    cutoffMethodDropdown.Value = 'Coefficient Threshold';
    cutoffMethodDropdown.Layout.Row = 1;
    cutoffMethodDropdown.Layout.Column = 2;
    cutoffMethodDropdown.ValueChangedFcn = @updateCutoffControls;
    
    % Cutoff value
    cutoffValueLabel = uilabel(cutoffGrid);
    cutoffValueLabel.Text = 'Threshold:';
    cutoffValueLabel.Layout.Row = 2;
    cutoffValueLabel.Layout.Column = 1;
    
    cutoffValueEdit = uieditfield(cutoffGrid, 'numeric');
    cutoffValueEdit.Value = 0.01;  % Default coefficient threshold
    cutoffValueEdit.Layout.Row = 2;
    cutoffValueEdit.Layout.Column = 2;
    
    % Run button
    runButton = uibutton(paramGrid);
    runButton.Text = 'Run LASSO Analysis';
    runButton.BackgroundColor=[.2 .2 .2];
    runButton.FontColor='w';
    runButton.Layout.Row = 9;
    runButton.Layout.Column = 1;
    runButton.ButtonPushedFcn = @runLassoAnalysis;
    
    % Results Panel Components - Fix for tab group layout
    % Create a grid layout inside the results panel to ensure proper expansion
    resultsGrid = uigridlayout(resultsPanel, [1 1]);
    resultsGrid.RowHeight = {'1x'};
    resultsGrid.ColumnWidth = {'1x'};
    
    % Create tab group inside the grid layout to ensure it fills the space
    resultsTabs = uitabgroup(resultsGrid);
    resultsTabs.Layout.Row = 1;
    resultsTabs.Layout.Column = 1;
    
    % Coefficient Path Tab
    coeffPathTab = uitab(resultsTabs);
    coeffPathTab.Title = 'Coefficient Path';
    coeffAxesGrid = uigridlayout(coeffPathTab, [1 1]);
    coeffAxesGrid.RowHeight = {'1x'};
    coeffAxesGrid.ColumnWidth = {'1x'};
    
    coeffAxes = uiaxes(coeffAxesGrid);
    coeffAxes.Layout.Row = 1;
    coeffAxes.Layout.Column = 1;
    coeffAxes.XLabel.String = 'Log(Lambda)';
    coeffAxes.YLabel.String = 'Coefficients';
    title(coeffAxes, 'LASSO Coefficient Path');
    
    % CV Plot Tab
    cvPlotTab = uitab(resultsTabs);
    cvPlotTab.Title = 'Cross-Validation';
    cvAxesGrid = uigridlayout(cvPlotTab, [1 1]);
    cvAxesGrid.RowHeight = {'1x'};
    cvAxesGrid.ColumnWidth = {'1x'};
    
    cvAxes = uiaxes(cvAxesGrid);
    cvAxes.Layout.Row = 1;
    cvAxes.Layout.Column = 1;
    cvAxes.XLabel.String = 'Log(Lambda)';
    cvAxes.YLabel.String = 'Mean Squared Error';
    title(cvAxes, 'Cross-Validation Error');
    
    % Coefficients Table Tab
    coeffTableTab = uitab(resultsTabs);
    coeffTableTab.Title = 'Coefficients';
    coeffTableGrid = uigridlayout(coeffTableTab, [1 1]);
    coeffTableGrid.RowHeight = {'1x'};
    coeffTableGrid.ColumnWidth = {'1x'};
    
    coeffTable = uitable(coeffTableGrid);
    coeffTable.Layout.Row = 1;
    coeffTable.Layout.Column = 1;
    coeffTable.ColumnName = {'Variable', 'Coefficient'};
    
    % Selected Variables Tab - NEW
    selectedVarsTab = uitab(resultsTabs);
    selectedVarsTab.Title = 'Selected Variables';
    selectedVarsGrid = uigridlayout(selectedVarsTab, [2 3]);
    selectedVarsGrid.RowHeight = {'1x', 24};
    selectedVarsGrid.ColumnWidth = {'1x','fit','1x'};
    
    selectedVarsTable = uitable(selectedVarsGrid);
    selectedVarsTable.Layout.Row = 1;
    selectedVarsTable.Layout.Column = [1 3];
    selectedVarsTable.ColumnName = {'Rank', 'Variable', 'Coefficient'};

    saveSelectedVarsBtn = uibutton(selectedVarsGrid);
    saveSelectedVarsBtn.Text = 'Save Selected Variables';
    saveSelectedVarsBtn.BackgroundColor=[.2 .2 .2];
    saveSelectedVarsBtn.FontColor='w';    
    saveSelectedVarsBtn.Layout.Row = 2;
    saveSelectedVarsBtn.Layout.Column = 2;
    saveSelectedVarsBtn.ButtonPushedFcn = @saveSelectedVars;
    
    
    % Save model tab
    saveModelTab = uitab(resultsTabs);
    saveModelTab.Title = 'Save Model';
    saveModelGrid = uigridlayout(saveModelTab, [3 3]);
    saveModelGrid.RowHeight = {30, 24, '1x'};
    saveModelGrid.ColumnWidth = {'fit', 'fit', '1x'};
    
    % Save model controls
    modelNameLabel = uilabel(saveModelGrid);
    modelNameLabel.Text = 'Model Name:';
    modelNameLabel.Layout.Row = 1;
    modelNameLabel.Layout.Column = 1;
    
    modelNameEdit = uieditfield(saveModelGrid, 'text');
    modelNameEdit.Value = 'LassoModel';
    modelNameEdit.Layout.Row = 1;
    modelNameEdit.Layout.Column = [2 3];
    
    saveModelBtn = uibutton(saveModelGrid);
    saveModelBtn.Text = 'Save Model';
    saveModelBtn.BackgroundColor=[.2 .2 .2];
    saveModelBtn.FontColor='w';    
        
    saveModelBtn.Layout.Row = 2;
    saveModelBtn.Layout.Column = 2;
    saveModelBtn.ButtonPushedFcn = @saveModel;
    
    % Store application data
    appData = struct();
    appData.X = X;
    appData.y = y;
    appData.varNames = varNames;
    appData.B = [];
    appData.FitInfo = [];
    appData.LassoModel = [];
    appData.SelectedVars = [];
    
    % Set app data
    fig.UserData = appData;
    drawnow;
    fig.Visible = 'on';
    
    % Callback Functions
    function updateLambdaControls(src, ~)
        switch src.Value
            case 'Manual'
                lambdaEdit.Enable = 'on';
                cvFoldsSpinner.Enable = 'off';
            case {'CV', 'MSE'}
                lambdaEdit.Enable = 'off';
                cvFoldsSpinner.Enable = 'on';
        end
    end
    
    function updateCutoffControls(src, ~)
        switch src.Value
            case 'Coefficient Threshold'
                cutoffValueLabel.Text = 'Threshold:';
                cutoffValueEdit.Value = 0.01;
            case 'Top N Variables'
                cutoffValueLabel.Text = 'Top N:';
                cutoffValueEdit.Value = min(10, n_vars);
        end
    end
    
    function runLassoAnalysis(~, ~)
        fw = gui.myWaitbar(fig);
        appData = fig.UserData;
        
        try
            % Get LASSO parameters
            alpha = alphaSlider.Value;
            numLambda = numLambdaSpinner.Value;
            standardize = standardizeCheck.Value;
            
            % Configure lambda selection method
            lambdaMethod = lambdaMethodDropdown.Value;
            
            switch lambdaMethod
                case 'CV'
                    kFolds = cvFoldsSpinner.Value;
                    [B, FitInfo] = lasso(X, y, 'Alpha', alpha, 'CV', kFolds, ...
                        'NumLambda', numLambda, 'Standardize', standardize);
                    lambdaOpt = FitInfo.LambdaMinMSE;
                    idxOpt = FitInfo.IndexMinMSE;
                    
                    % Plot the CV curve
                    plot(cvAxes, log(FitInfo.Lambda), FitInfo.MSE);
                    hold(cvAxes, 'on');
                    plot(cvAxes, log(lambdaOpt)*[1,1], ylim(cvAxes), 'r--');
                    hold(cvAxes, 'off');
                    title(cvAxes, 'Cross-Validation Mean Squared Error');
                    xlabel(cvAxes, 'Log(Lambda)');
                    ylabel(cvAxes, 'Mean Squared Error');
                    grid(cvAxes, 'on');
                    
                case 'MSE'
                    [B, FitInfo] = lasso(X, y, 'Alpha', alpha, ...
                        'NumLambda', numLambda, 'Standardize', standardize);
                    
                    % Choose lambda that gives minimum MSE on training data
                    yhat = X * B + FitInfo.Intercept(ones(1,length(FitInfo.Lambda)));
                    mse = mean((repmat(y, 1, length(FitInfo.Lambda)) - yhat).^2);
                    [~, idxOpt] = min(mse);
                    lambdaOpt = FitInfo.Lambda(idxOpt);
                    
                case 'Manual'
                    lambdaOpt = lambdaEdit.Value;
                    
                    % Generate a sequence of lambda values ending with the user's choice
                    lambdaValues = exp(linspace(log(lambdaOpt*100), log(lambdaOpt), numLambda));
                    
                    [B, FitInfo] = lasso(X, y, 'Alpha', alpha, ...
                        'Lambda', lambdaValues, 'Standardize', standardize);
                    
                    idxOpt = numLambda; % Last lambda is the one specified by user
            end
            
            % Plot coefficient paths
            plot(coeffAxes, log(FitInfo.Lambda), B);
            hold(coeffAxes, 'on');
            if exist('lambdaOpt', 'var')
                plot(coeffAxes, log(lambdaOpt)*[1,1], ylim(coeffAxes), 'r--');
            end
            hold(coeffAxes, 'off');
            title(coeffAxes, 'LASSO Coefficient Paths');
            xlabel(coeffAxes, 'Log(Lambda)');
            ylabel(coeffAxes, 'Coefficients');
            grid(coeffAxes, 'on');
            
            % Create a legend for coefficient paths
            if length(varNames) < 10
                legend(coeffAxes, varNames, 'Location', 'eastoutside');
            else
                legend(coeffAxes, varNames(1:10), 'Location', 'eastoutside');
            end
            
            % Update coefficients table
            if exist('idxOpt', 'var')
                coefs = [FitInfo.Intercept(idxOpt); B(:, idxOpt)];
                coefNames = ['(Intercept)'; varNames(:)];
                nonZeroIdx = abs(coefs) > 1e-10;
                tableData = table(coefNames(nonZeroIdx), coefs(nonZeroIdx), ...
                    'VariableNames', {'Variable', 'Coefficient'});
                
                % Sort by absolute coefficient value
                [~, sortIdx] = sort(abs(tableData.Coefficient), 'descend');
                tableData = tableData(sortIdx, :);
                assignin('base',"tableData",tableData)
                tableData.Variable = cellstr(tableData.Variable);
                coeffTable.Data = table2cell(tableData);
                coeffTable.ColumnName = {'Variable', 'Coefficient'};
                
                % Apply variable selection cutoff
                cutoffMethod = cutoffMethodDropdown.Value;
                cutoffValue = cutoffValueEdit.Value;
                
                % Extract model coefficients (without intercept)
                modelCoefs = B(:, idxOpt);
                absCoefs = abs(modelCoefs);
                
                % Apply cutoff based on the selected method
                switch cutoffMethod
                    case 'Coefficient Threshold'
                        selectedIdx = absCoefs >= cutoffValue;
                    case 'Top N Variables'
                        % Ensure N is within bounds
                        N = min(max(1, round(cutoffValue)), sum(absCoefs > 0));
                        [~, sortedIdx] = sort(absCoefs, 'descend');
                        selectedIdx = false(size(modelCoefs));
                        selectedIdx(sortedIdx(1:N)) = true;
                end
                
                % Create table with selected variables
                if any(selectedIdx)
                    selCoefs = modelCoefs(selectedIdx);
                    selVarNames = varNames(selectedIdx);
                    [~, rankIdx] = sort(abs(selCoefs), 'descend');
                    
                    %assignin('base',"selCoefs",selCoefs)
                    %assignin('base',"selVarNames",selVarNames)

                    selTable = table((1:length(selCoefs))', ...
                                     selVarNames(rankIdx), ...
                                     selCoefs(rankIdx), ...
                                     'VariableNames', {'Rank', 'Variable', 'Coefficient'});
                    
                    %assignin('base', "selTable", selTable);

                    selTable.Variable = cellstr(selTable.Variable);                    
                    selectedVarsTable.Data = table2cell(selTable);
                    selectedVarsTable.ColumnName = {'Rank', 'Variable', 'Coefficient'};
                    
                    % Store selected variables
                    appData.SelectedVars = selVarNames(rankIdx);
                else
                    % No variables selected
                    selectedVarsTable.Data = {};
                    appData.SelectedVars = {};
                    
                    uialert(fig, 'No variables meet the selection criteria. Consider adjusting the threshold.', ...
                           'No Variables Selected', 'Icon', 'warning');
                end
            end
            
            % Store the results in app data
            appData.B = B;
            appData.FitInfo = FitInfo;
            if exist('idxOpt', 'var')
                appData.LassoModel = struct('Intercept', FitInfo.Intercept(idxOpt), ...
                    'Coefficients', B(:, idxOpt), ...
                    'Lambda', FitInfo.Lambda(idxOpt), ...
                    'Alpha', alpha, ...
                    'PredictorVars', {varNames}, ...
                    'Standardize', standardize);
            end
            fig.UserData = appData;
            gui.myWaitbar(fig, fw);    
        catch ME
            gui.myWaitbar(fig, fw, true);
            uialert(fig, ['Error in LASSO analysis: ', ME.message], 'Analysis Error', 'Icon', 'error');
        end
        
    end
    
    function saveSelectedVars(~, ~)
        appData = fig.UserData;
        
        % Check if we have a model to save
        if ~isfield(appData, 'LassoModel') || isempty(appData.LassoModel)
            uialert(fig, 'Please run the LASSO analysis first.', 'No Model', 'Icon', 'warning');
            return;
        end
        
        try
            % Get the model name
            modelName = modelNameEdit.Value;
            if isempty(modelName)
                modelName = 'SelectedGenes';
            end
            
            % Create a variable with the model name in the workspace
            model = appData.LassoModel;
            
            % Add selected variables to the model
            model.SelectedVariables = appData.SelectedVars;
            
            assignin('base', "SelectedGenes", model.SelectedVariables);
            
            % Optionally save to file
            %[file, path] = uiputfile({'*.mat', 'MATLAB Data File (*.mat)'}, ...
            %    'Save LASSO Model', [modelName, '.mat']);
            %if ~isequal(file, 0)
                %fullPath = fullfile(path, file);
                %save(fullPath, 'model');
                %uialert(fig, ['Model saved to workspace as "', modelName, '" and to file "', fullPath, '".'], ...
                %    'Model Saved', 'Icon', 'success');
            %else
                uialert(fig, ['Selected genes saved to workspace as "', modelName, '".'], ...
                    'Genes Saved', 'Icon', 'success');
            %end
            
        catch ME
            uialert(fig, ['Error saving model: ', ME.message], 'Save Error', 'Icon', 'error');
        end
    end


    function saveModel(~, ~)
        appData = fig.UserData;
        
        % Check if we have a model to save
        if ~isfield(appData, 'LassoModel') || isempty(appData.LassoModel)
            uialert(fig, 'Please run the LASSO analysis first.', 'No Model', 'Icon', 'warning');
            return;
        end
        
        try
            % Get the model name
            modelName = modelNameEdit.Value;
            if isempty(modelName)
                modelName = 'LassoModel';
            end
            
            % Create a variable with the model name in the workspace
            model = appData.LassoModel;
            
            % Add selected variables to the model
            model.SelectedVariables = appData.SelectedVars;
            
            assignin('base', modelName, model);
            
            % Optionally save to file
            %[file, path] = uiputfile({'*.mat', 'MATLAB Data File (*.mat)'}, ...
            %    'Save LASSO Model', [modelName, '.mat']);
            %if ~isequal(file, 0)
                %fullPath = fullfile(path, file);
                %save(fullPath, 'model');
                %uialert(fig, ['Model saved to workspace as "', modelName, '" and to file "', fullPath, '".'], ...
                %    'Model Saved', 'Icon', 'success');
            %else
                uialert(fig, ['Model saved to workspace as "', modelName, '".'], ...
                    'Model Saved', 'Icon', 'success');
            %end
            
        catch ME
            uialert(fig, ['Error saving model: ', ME.message], 'Save Error', 'Icon', 'error');
        end
    end

end
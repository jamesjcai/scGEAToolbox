classdef DEResultViewApp_simple < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure              matlab.ui.Figure
        TabGroup              matlab.ui.container.TabGroup
        FullTableTab          matlab.ui.container.Tab
        UpRegulatedTab        matlab.ui.container.Tab
        DownRegulatedTab      matlab.ui.container.Tab
        SignificanceTab       matlab.ui.container.Tab
        VolcanoPlotTab        matlab.ui.container.Tab
        
        % Tables
        FullTable             matlab.ui.control.Table
        UpRegulatedTable      matlab.ui.control.Table
        DownRegulatedTable    matlab.ui.control.Table
        
        % Significance controls
        PValueSlider          matlab.ui.control.Slider
        PValueLabel           matlab.ui.control.Label
        PValueEditField       matlab.ui.control.NumericEditField
        FCSlider              matlab.ui.control.Slider
        FCLabel               matlab.ui.control.Label
        FCEditField           matlab.ui.control.NumericEditField
        ApplyButton           matlab.ui.control.Button
        SignificantGenesLabel matlab.ui.control.Label
        
        % Volcano plot
        VolcanoPanel          matlab.ui.container.Panel
        VolcanoAxis           matlab.ui.control.UIAxes
    end
    
    % App data
    properties (Access = private)
        DEData % Differential expression data
        PValueColumn = 'pvalue'; % Default column names - can be changed
        LogFCColumn = 'logFC';
        GeneColumn = 'GeneID';
        
        % Threshold values
        PValueThreshold = 0.05;
        FCThreshold = 1.0; % logFC threshold
    end
    
    methods (Access = private)
        
        function updateTables(app)
            % Update the full table display
            app.FullTable.Data = app.DEData;
            
            % Extract up-regulated genes based on thresholds
            upRegIndices = app.DEData.(app.LogFCColumn) >= app.FCThreshold & ...
                          app.DEData.(app.PValueColumn) <= app.PValueThreshold;
            upRegData = app.DEData(upRegIndices, :);
            app.UpRegulatedTable.Data = upRegData;
            
            % Extract down-regulated genes based on thresholds
            downRegIndices = app.DEData.(app.LogFCColumn) <= -app.FCThreshold & ...
                            app.DEData.(app.PValueColumn) <= app.PValueThreshold;
            downRegData = app.DEData(downRegIndices, :);
            app.DownRegulatedTable.Data = downRegData;
            
            % Update the significant genes count label
            app.SignificantGenesLabel.Text = sprintf('Significant Genes: %d (Up: %d, Down: %d)', ...
                sum(upRegIndices) + sum(downRegIndices), sum(upRegIndices), sum(downRegIndices));
        end
        
        function updateVolcanoPlot(app)
            % Clear the axis
            cla(app.VolcanoAxis);
            
            % Get the p-values and logFC values
            pvals = app.DEData.(app.PValueColumn);
            logfc = app.DEData.(app.LogFCColumn);
            
            % Calculate -log10(p-value) for the y-axis
            negLogPvals = -log10(pvals);
            
            % Create scatter categorization
            upRegIndices = logfc >= app.FCThreshold & pvals <= app.PValueThreshold;
            downRegIndices = logfc <= -app.FCThreshold & pvals <= app.PValueThreshold;
            nonsigIndices = ~(upRegIndices | downRegIndices);
            
            % Plot non-significant points (gray)
            scatter(app.VolcanoAxis, logfc(nonsigIndices), negLogPvals(nonsigIndices), 20, [0.7 0.7 0.7], 'filled');
            hold(app.VolcanoAxis, 'on');
            
            % Plot up-regulated points (red)
            scatter(app.VolcanoAxis, logfc(upRegIndices), negLogPvals(upRegIndices), 30, 'r', 'filled');
            
            % Plot down-regulated points (blue)
            scatter(app.VolcanoAxis, logfc(downRegIndices), negLogPvals(downRegIndices), 30, 'b', 'filled');
            
            % Add threshold lines
            plot(app.VolcanoAxis, [app.FCThreshold app.FCThreshold], ylim(app.VolcanoAxis), 'k--');
            plot(app.VolcanoAxis, [-app.FCThreshold -app.FCThreshold], ylim(app.VolcanoAxis), 'k--');
            plot(app.VolcanoAxis, xlim(app.VolcanoAxis), [-log10(app.PValueThreshold) -log10(app.PValueThreshold)], 'k--');
            
            % Add labels and legend
            xlabel(app.VolcanoAxis, 'log2 Fold Change');
            ylabel(app.VolcanoAxis, '-log10(p-value)');
            title(app.VolcanoAxis, 'Volcano Plot');
            legend(app.VolcanoAxis, {'Non-significant', 'Up-regulated', 'Down-regulated'}, 'Location', 'best');
            
            hold(app.VolcanoAxis, 'off');
            
            % Add text annotations for the top significant genes
            numToLabel = min(10, sum(upRegIndices | downRegIndices));
            if numToLabel > 0 && isfield(app.DEData, app.GeneColumn)
                % Combine up and down regulated
                sigIndices = upRegIndices | downRegIndices;
                sigPvals = pvals(sigIndices);
                sigLogFC = logfc(sigIndices);
                sigNegLogPvals = negLogPvals(sigIndices);
                sigGenes = app.DEData.(app.GeneColumn)(sigIndices);
                
                % Sort by significance
                [~, sortIdx] = sort(sigPvals);
                topIndices = sortIdx(1:numToLabel);
                
                % Add text labels
                hold(app.VolcanoAxis, 'on');
                for i = 1:numToLabel
                    idx = topIndices(i);
                    text(app.VolcanoAxis, sigLogFC(idx), sigNegLogPvals(idx), ...
                         char(sigGenes(idx)), 'FontSize', 8, 'VerticalAlignment', 'bottom');
                end
                hold(app.VolcanoAxis, 'off');
            end
        end
    end
    
    % Callbacks
    methods (Access = private)
        
        function pValueSliderChanged(app, ~)
            % Update the p-value threshold
            app.PValueThreshold = app.PValueSlider.Value;
            app.PValueEditField.Value = app.PValueThreshold;
            
            % Update tables and volcano plot
            updateTables(app);
            updateVolcanoPlot(app);
        end
        
        function pValueEditFieldChanged(app, ~)
            % Update the p-value threshold
            app.PValueThreshold = app.PValueEditField.Value;
            app.PValueSlider.Value = app.PValueThreshold;
            
            % Update tables and volcano plot
            updateTables(app);
            updateVolcanoPlot(app);
        end
        
        function fcSliderChanged(app, ~)
            % Update the fold change threshold
            app.FCThreshold = app.FCSlider.Value;
            app.FCEditField.Value = app.FCThreshold;
            
            % Update tables and volcano plot
            updateTables(app);
            updateVolcanoPlot(app);
        end
        
        function fcEditFieldChanged(app, ~)
            % Update the fold change threshold
            app.FCThreshold = app.FCEditField.Value;
            app.FCSlider.Value = app.FCThreshold;
            
            % Update tables and volcano plot
            updateTables(app);
            updateVolcanoPlot(app);
        end
        
        function applyButtonPushed(app, ~)
            % Apply all threshold changes
            updateTables(app);
            updateVolcanoPlot(app);
        end
    end
    
    % App creation and setup
    methods (Access = private)
        
        function createComponents(app)
            % Create the figure and components
            app.UIFigure = uifigure('Name', 'DE Result View', 'Position', [100, 100, 800, 600]);
            
            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure, 'Position', [1, 1, 800, 600]);
            
            % Full Table Tab
            app.FullTableTab = uitab(app.TabGroup, 'Title', 'Full DE Table');
            app.FullTable = uitable(app.FullTableTab, 'Position', [20, 20, 760, 540]);
            app.FullTable.ColumnSortable = true;
            
            % Up-Regulated Tab
            app.UpRegulatedTab = uitab(app.TabGroup, 'Title', 'Up-Regulated Genes');
            app.UpRegulatedTable = uitable(app.UpRegulatedTab, 'Position', [20, 20, 760, 540]);
            app.UpRegulatedTable.ColumnSortable = true;
            
            % Down-Regulated Tab
            app.DownRegulatedTab = uitab(app.TabGroup, 'Title', 'Down-Regulated Genes');
            app.DownRegulatedTable = uitable(app.DownRegulatedTab, 'Position', [20, 20, 760, 540]);
            app.DownRegulatedTable.ColumnSortable = true;
            
            % Significance Threshold Tab
            app.SignificanceTab = uitab(app.TabGroup, 'Title', 'Significance Thresholds');
            
            % P-value controls
            app.PValueLabel = uilabel(app.SignificanceTab, 'Position', [50, 500, 150, 22], 'Text', 'P-value threshold:');
            app.PValueSlider = uislider(app.SignificanceTab, 'Position', [200, 510, 400, 3]);
            app.PValueSlider.Limits = [0.001, 0.1];
            app.PValueSlider.Value = app.PValueThreshold;
            app.PValueSlider.ValueChangedFcn = @(src, event) pValueSliderChanged(app, event);
            
            app.PValueEditField = uieditfield(app.SignificanceTab, 'numeric', 'Position', [620, 500, 100, 22]);
            app.PValueEditField.Limits = [0.001, 0.1];
            app.PValueEditField.Value = app.PValueThreshold;
            app.PValueEditField.ValueChangedFcn = @(src, event) pValueEditFieldChanged(app, event);
            
            % Fold change controls
            app.FCLabel = uilabel(app.SignificanceTab, 'Position', [50, 450, 150, 22], 'Text', 'log2FC threshold:');
            app.FCSlider = uislider(app.SignificanceTab, 'Position', [200, 460, 400, 3]);
            app.FCSlider.Limits = [0, 4];
            app.FCSlider.Value = app.FCThreshold;
            app.FCSlider.ValueChangedFcn = @(src, event) fcSliderChanged(app, event);
            
            app.FCEditField = uieditfield(app.SignificanceTab, 'numeric', 'Position', [620, 450, 100, 22]);
            app.FCEditField.Limits = [0, 4];
            app.FCEditField.Value = app.FCThreshold;
            app.FCEditField.ValueChangedFcn = @(src, event) fcEditFieldChanged(app, event);
            
            % Apply button
            app.ApplyButton = uibutton(app.SignificanceTab, 'Position', [350, 380, 100, 30], 'Text', 'Apply');
            app.ApplyButton.ButtonPushedFcn = @(src, event) applyButtonPushed(app, event);
            
            % Significant genes count label
            app.SignificantGenesLabel = uilabel(app.SignificanceTab, 'Position', [250, 330, 300, 22], 'Text', 'Significant Genes: 0 (Up: 0, Down: 0)');
            app.SignificantGenesLabel.HorizontalAlignment = 'center';
            app.SignificantGenesLabel.FontWeight = 'bold';
            
            % Volcano Plot Tab
            app.VolcanoPlotTab = uitab(app.TabGroup, 'Title', 'Volcano Plot');
            app.VolcanoPanel = uipanel(app.VolcanoPlotTab, 'Position', [20, 20, 760, 540]);
            app.VolcanoAxis = uiaxes(app.VolcanoPanel, 'Position', [30, 30, 700, 480]);
        end
    end
    
    % Public methods
    methods (Access = public)
        
        function loadData(app, data)
            % Load differential expression data into the app
            app.DEData = data;
            
            % Try to auto-detect column names
            colNames = data.Properties.VariableNames;
            
            % Look for p-value columns
            pvalCols = contains(lower(colNames), {'pval', 'padj', 'fdr', 'qval','p_val_adj'});
            if any(pvalCols)
                app.PValueColumn = colNames{find(pvalCols, 1)};
            end
            
            % Look for logFC columns
            logFCCols = contains(lower(colNames), {'logfc', 'log2fc', 'log_fc', 'fold','avg_log2FC'});
            if any(logFCCols)
                app.LogFCColumn = colNames{find(logFCCols, 1)};
            end
            
            % Look for gene ID columns
            geneCols = contains(lower(colNames), {'gene', 'id', 'symbol', 'name'});
            if any(geneCols)
                app.GeneColumn = colNames{find(geneCols, 1)};
            end
            
            % Set up tables
            app.FullTable.ColumnName = colNames;
            app.UpRegulatedTable.ColumnName = colNames;
            app.DownRegulatedTable.ColumnName = colNames;
            
            % Update display
            updateTables(app);
            updateVolcanoPlot(app);
        end
    end
    
    % App initialization and construction
    methods (Access = public)
        
        function app = DEResultViewApp_simple(data)
            % Create the components
            createComponents(app);
            
            % Load data if provided
            if nargin > 0 && ~isempty(data)
                loadData(app, data);
            end
            
            % Show the figure
            app.UIFigure.Visible = 'on';
        end
    end
end

% Example usage of the DEResultViewApp
function sampleUsage()
    % Create a sample DE results table
    numGenes = 1000;
    genes = cellstr(strcat('Gene', num2str((1:numGenes)')));
    logFC = randn(numGenes, 1) * 2;  % Random fold changes
    pvalues = rand(numGenes, 1) .^ 2; % Random p-values (skewed toward lower values)
    
    % Create a table
    deData = table(genes, logFC, pvalues, 'VariableNames', {'GeneID', 'logFC', 'pvalue'});
    
    % Create and display the app
    app = DEResultViewApp(deData);
end
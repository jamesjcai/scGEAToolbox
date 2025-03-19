function sc_pseudotimeexplorer(X, genelist, s, varargin)
%load cleandata.mat
%s=s_tsne;
parentfig = []; 
p = inputParser;
addRequired(p, 'X', @isnumeric);
addRequired(p, 'genelist', @isstring);
addRequired(p, 's', @isnumeric);
addOptional(p, 'method', "splinefit", @(x) (isstring(x) | ischar(x)) ...
    & ismember(lower(string(x)), "splinefit"));
addOptional(p, 'dim', 1, @isnumeric);
parse(p, X, genelist, s, varargin{:});
%method=p.Results.method;
dim = p.Results.dim;


%global psexplorer_timeid
%psexplorer_timeid=1;
hx=gui.myFigure;
hFig=hx.FigHandle;
hAx = axes('Parent', hFig);

if size(s, 2) >= 3
    hs = scatter3(hAx, s(:, 1), s(:, 2), s(:, 3), 10);
elseif size(s, 2) == 2
    hs = scatter(hAx, s(:, 1), s(:, 2), 10);
end

hx.addCustomButton('off',  @drawtrajectory, 'plotpicker-plot.gif', 'Plot pseudotime trajectory');
hx.addCustomButton('off',  @deleteselectedcells, 'plotpicker-scatter.gif', 'Delet selected cells');

hx.show();


    % function selectdimension(~, ~)
    %     dim = dim + 1;
    %     if dim > 3, dim = 1; end
    %     msgbox(sprintf('dim=%d\n', dim), ...
    %         'Set Dimension');
    % end


        function deleteselectedcells(~, ~)
            data = hs.BrushData;
            ptsSelected = find(data);
            X(:, ptsSelected) = [];
            s(ptsSelected, :) = [];
            [a, b] = view();
            if size(s, 2) >= 3
                hs = scatter3(hAx, s(:, 1), s(:, 2), s(:, 3), 10);
            elseif size(s, 2) == 2
                hs = scatter(hAx, s(:, 1), s(:, 2), 10);
            end
            view(a, b);
        end


        function drawtrajectory(~, ~)
            [t, xyz1] = pkg.i_pseudotime_by_splinefit(s, dim, false);
            hold on
            if size(xyz1, 2) >= 3
                plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-r', 'linewidth', 2);
                text(xyz1(1, 1), xyz1(1, 2), xyz1(1, 3), 'Start', ...
                    'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
                text(xyz1(end, 1), xyz1(end, 2), xyz1(end, 3), 'End', ...
                    'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            elseif size(xyz1, 2) == 2
                plot(xyz1(:, 1), xyz1(:, 2), '-r', 'linewidth', 2);
                text(xyz1(1, 1), xyz1(1, 2), 'Start', ...
                    'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
                text(xyz1(end, 1), xyz1(end, 2), 'End', ...
                    'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
            end
            hold off

            labels = {'Save expression X to variable named:', ...
                'Save pseudotime T to variable named:', ...
                'Save embedding S to variable named:'};
            vars = {'X_psexplorer', 't_psexplorer', 's_psexplorer'};
            values = {X, t, s};
            if ~(ismcc || isdeployed)
                msgfig = export2wsdlg(labels, vars, values);
                %         assignin('base',sprintf('psexplorerT%d',...
                %                  psexplorer_timeid),t);
                uiwait(msgfig)
            end
            answer = gui.myQuestdlg(parentfig, 'View expression of selected genes', ...
                'Pseudotime Function', ...
                'Yes', 'No', 'Yes');
            if isempty(answer), return; end
            switch answer
                case 'Yes'
                    r = corr(t, X', 'type', 'spearman'); % Calculate linear correlation between gene expression profile and T
                    [~, idxp] = maxk(r, 4); % Select top 4 positively correlated genes
                    [~, idxn] = mink(r, 3); % Select top 3 negatively correlated genes
                    selectedg = genelist([idxp, idxn]);
                    figure;
                    pkg.i_plot_pseudotimeseries(log2(X+1), ...
                        genelist, t, selectedg);
                case 'No'
                    return;
            end
        end

end


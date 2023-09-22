function sc_markerexplorer(X, genelist, s, varargin)
%load cleandata.mat
%s=s_tsne;

p = inputParser;
addRequired(p, 'X', @isnumeric);
addRequired(p, 'genelist', @isstring);
addRequired(p, 's', @isnumeric);
addOptional(p, 'method', "ranksum", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["ranksum", "mast"]));
addOptional(p, 'numfig', 1, @isnumeric);
parse(p, X, genelist, s, varargin{:});
method = p.Results.method;
numfig = p.Results.numfig;

mfolder = fileparts(mfilename('fullpath'));

%global mkexplorer_clustid
mkexplorer_clustid = 0;
hFig = figure('Name', 'Marker Gene Explorer');
hAx = axes('Parent', hFig);

if size(s, 2) >= 3
    scatter3(hAx, s(:, 1), s(:, 2), s(:, 3), 10);
elseif size(s, 2) == 2
    scatter(hAx, s(:, 1), s(:, 2), 10);
end

tb = uitoolbar(hFig);
tt = uitoggletool(tb, 'Separator', 'on');
% [img,map] = imread(fullfile(matlabroot,...
%             'toolbox','matlab','icons','tool_ellipse.gif'));
[img, map] = imread(fullfile(mfolder, ...
    'resources', 'brush.gif'));
ptImage = ind2rgb(img, map);
tt.CData = ptImage;
tt.Tooltip = 'Click and then brush/select cells';
tt.ClickedCallback = @BrushSwitcher;

    function BrushSwitcher(src, ~)
        state = src.State;
        if strcmp(state, 'on')
            hBr.Enable = 'on';
            %tt.CData = zeros(16,16,3);
        else
            hBr.Enable = 'off';
            %tt.CData = ptImage;
        end
end


    pt = uipushtool(tb, 'Separator', 'on');
    % [img,map] = imread(fullfile(matlabroot,...
    %             'toolbox','matlab','icons','profiler.gif'));
    [img, map] = imread(fullfile(mfolder, ...
        'resources', 'list.gif'));
    ptImage = ind2rgb(img, map);

    % defaultToolbar = findall(hFig,'Type','uitoolbar');
    % pt = uipushtool(defaultToolbar);
    % ptImage = rand(16,16,3);
    pt.CData = ptImage;
    pt.Tooltip = 'Select a gene to show expression';
    pt.ClickedCallback = @showmkgene;

        function showmkgene(~, ~)
            gsorted = sort(genelist);
            [indx, tf] = listdlg('PromptString', {'Select a gene', ...
                '', ''}, 'SelectionMode', 'single', 'ListString', gsorted);
            if tf == 1
                figure;
                sc_scattermarker(X, genelist, gsorted(indx), s, 5);
                %       [ax,bx]=view(); sc_markerscatter(X,genelist,gsorted(indx),s,3);
                %       view(ax,bx);
            end
    end


        pt3 = uipushtool(tb, 'Separator', 'on');
        % [img,map] = imread(fullfile(matlabroot,...
        %             'toolbox','matlab','icons','plotpicker-scatter.gif'));
        [img, map] = imread(fullfile(matlabroot, ...
            'toolbox', 'matlab', 'icons', 'plotpicker-stairs.gif'));

        ptImage = ind2rgb(img, map);
        pt3.CData = ptImage;
        pt3.Tooltip = 'Set numfig';
        pt3.ClickedCallback = @selectnumfig;

            function selectnumfig(~, ~)
                numfig = numfig + 1;
                if numfig > 10, numfig = 1; end
                msgbox(sprintf('numfig=%d\n', numfig), ...
                    'Set Number of Figures');
        end


            %h = uicontrol('Position',[5 5 150 30],'String','Calculate xdiff',...
            %              'Callback', @JCal);

            hBr = brush(hFig);
            % hBr.Enable='on';
            % hBr.Color = 'green';
            %hBr.ActionPostCallback = {@onBrushAction,X,genelist,s,method,numfig};
            hBr.ActionPostCallback = @onBrushAction;

            gui.add_3dcamera(tb);


            % ref: https://www.mathworks.com/matlabcentral/answers/385226-how-to-use-the-data-brush-tool-to-automatically-save-selected-points-in-multiple-line-plots

                function onBrushAction(~, eventdata)
                    % global mkexplorer_clustid
                    % Extract plotted graphics objects
                    % Invert order because "Children" property is in reversed plotting order
                    hLines = flipud(eventdata.Axes.Children);

                    % Loop through each graphics object
                    for k = 1:1 %numel(hLines)
                        % Check that the property is valid for that type of object
                        % Also check if any points in that object are selected
                        if isprop(hLines(k), 'BrushData') && any(hLines(k).BrushData)
                            % Output the selected data to the base workspace with assigned name
                            ptsSelected = logical(hLines(k).BrushData.');
                            switch method
                                case 'mast'
                                    [markerlist] = sc_pickmarkers2(X, genelist, 1+ptsSelected, 2);
                                case 'ranksum'
                                    [markerlist] = sc_pickmarkers(X, genelist, 1+ptsSelected, 2);
                            end
                            [ax, bx] = view();
                            for kkk = 1:numfig
                                figure;
                                for kk = 1:min([9, length(markerlist)])
                                    subplot(3, 3, kk)
                                    sc_scattermarker(X, genelist, s, ...
                                        markerlist(kk+9*(kkk - 1)), 3, 3, false);
                                    view(ax, bx);
                                end
                            end
                            mkexplorer_clustid = mkexplorer_clustid + 1;
                            assignin('base', sprintf('mkexplorerL%d', ...
                                mkexplorer_clustid), markerlist);
                        end
                    end
            end

            end

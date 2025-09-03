function sc_uitabgrpfig_expplot(y, glist, s, parentfig, cazcel)

if nargin < 5, cazcel = []; end
if nargin < 4, parentfig = []; end
% if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
%     p = parentfig.Position;
%     cx = [p(1)+p(3)/2 p(2)+p(4)/2];
% end

% https://www.mathworks.com/help/rptgen/ug/compile-a-presentation-program.html
if ismcc || isdeployed, makePPTCompilable(); end
import mlreportgen.ppt.*;

% pw1 = fileparts(mfilename('fullpath'));
% pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');


hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
hFig.Position(3) = hFig.Position(3) * 1.8;

n = length(glist);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();
tab = cell(n,1);
% ax0 = cell(n,1);
ax = cell(n,2);

idx = 1;
focalg = glist(idx);

% y = cell(n,1);
% for k=1:n
%     y{k} = sce.X(sce.g == glist(k), :);
% end

for k = 1:n
    c = y{k};
    if issparse(c), c = full(c); end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s', glist(k)));
    
    %{
    t = tiledlayout(1,2,'Parent',tab{k});
    ax1 = nexttile;
    hpl{k,1} = scatter3(s(:,1), s(:,2), s(:,3), 5, c, 'filled','Parent', ax1);
    ax2 = nexttile;
    hpl{k,2} = scatter(s(:,1), s(:,2), 5, c, 'filled','Parent', ax2);
    %}
    
    % ax0{k} = axes('parent',tab{k});

    ax{k, 1} = subplot(1, 2, 1,'Parent', tab{k});
    % ax{k,1}.Tag = sprintf('axes_tab%d_left',k);

    if size(s,2)>2
        scatter3(ax{k,1}, s(:,1), s(:,2), s(:,3), 5, c, 'filled');
    else
        scatter(ax{k,1}, s(:,1), s(:,2), 5, c, 'filled');
    end
    if ~isempty(cazcel)
        view(ax{k,1}, [cazcel(1), cazcel(2)]);
    end
        title(ax{k,1}, glist(k));
        subtitle(ax{k,1}, gui.i_getsubtitle(c));

    gui.i_setautumncolor(c, a, true, any(c==0), ax{k,1}, parentfig);

    % ax{k,2} = subplot(1,2,2);

    ax{k,2} = subplot(1, 2, 2,'Parent', tab{k});
    % ax{k,2}.Tag = sprintf('axes_tab%d_right', k);   % 👈 assign unique Tag
    

        scatter(ax{k,2}, s(:,1), s(:,2), 5, c, 'filled');
        stem3(ax{k,2}, s(:,1), s(:,2), c, 'marker', 'none', 'color', 'm');
        hold(ax{k,2}, "on");
        scatter3(ax{k,2}, s(:,1), s(:,2), zeros(size(s(:,2))), 5, c, 'filled');
        
        title(ax{k,2}, glist(k));
        subtitle(ax{k,2}, gui.i_getsubtitle(c));

        axis(ax{k,2},'vis3d'); grid(ax{k,2},'on');

        % Enable rotate3d for the whole figure
        hRotate = rotate3d(hFig);
        set(hRotate,'Enable','on');

        % disp('rotate3d enabled');
        % Register callbacks
        % hRotate.ActionPreCallback  = @(src,evnt) startDrag(hFig, evnt, ax{k,2});
        
        %{
        ax = axes('Parent', tab{k});
        hold(ax, 'on');
        ax1 = subplot(1,2,1,tab{k});
        hold(ax1,'on')
        hpl{k,1} = scatter3(s(:,1), s(:,2), s(:,3), 5, c, 'filled','Parent', ax1);
        ax2 = subplot(1,2,2,tab{k});
        hold(ax2,'on')
        hpl{k,2} = scatter(s(:,1), s(:,2), 5, c, 'filled','Parent', ax2);
        %}
    
        %{
        hax{k} = axes('Parent', tab{k});
        if size(s,2)>=3
            hpl{k} = scatter3(s(:,1), s(:,2), s(:,3), 5, c, 'filled','Parent', hax{k});
        else
            hpl{k} = scatter(s(:,1), s(:,2), 5, c, 'filled','Parent', hax{k});
        end
        title(hax{k}, targetg(k));
        subtitle(hax{k}, gui.i_getsubtitle(c));
        %}
    
end
trackedAxes = [ax{:,1}, ax{:,2}];   % only the right ones, or ax(:) if you want all
hRotate.ActionPostCallback = @(src,evnt) stopDrag(evnt, trackedAxes);
% hRotate.ActionPostCallback = @(src,evnt) stopDrag(hFig, evnt);
  
tabgp.SelectionChangedFcn=@displaySelection;
hx.addCustomButton('off',  @in_genecards, 'www.jpg', 'GeneCards...');
hx.addCustomButton('off', @in_savedata, "floppy-disk-arrow-in.jpg", 'Save Gene List...');
hx.addCustomButton('off', {@gui.callback_RunGeneAgent, glist}, "mw-microprocessor.jpg", 'Run GeneAgent...');
hx.addCustomButton('off',  @in_mergetabs, 'Brightness-3--Streamline-Core.jpg', 'Show on the same figure...');

hx.show(parentfig)

    function in_mergetabs(~, ~)
        figure;
        for kx = 1:n
            hAx2 = nexttile;
            hAx1 = ax{kx,1};
            cloneAxes(hAx1, hAx2);
        end
        figure;
        for kx = 1:n
            hAx2 = nexttile;
            hAx1 = ax{kx,2};
            cloneAxes(hAx1, hAx2);
        end
    end

    function cloneAxes(hAx1, hAx2)
        copyobj(allchild(hAx1), hAx2);
        
        props = {'XLim','YLim','ZLim','XScale','YScale','ZScale',...
                 'XDir','YDir','ZDir','Colormap','CLim','View'};
        for kk = 1:numel(props)
            try
                set(hAx2, props{kk}, get(hAx1, props{kk}));
            catch
            end
        end
        
        xlabel(hAx2, get(get(hAx1,'XLabel'),'String'));
        ylabel(hAx2, get(get(hAx1,'YLabel'),'String'));
        title(hAx2,  get(get(hAx1,'Title'),'String'));
        subtitle(hAx2,  get(get(hAx1,'Subtitle'),'String'));

            gridProps = {'XGrid','YGrid','ZGrid', ...
             'XMinorGrid','YMinorGrid','ZMinorGrid', ...
             'Box'};
        for kk = 1:numel(gridProps)
            set(hAx2, gridProps{kk}, get(hAx1, gridProps{kk}));
        end

    end
    
    function in_savedata(~,~)
        gui.i_exporttable(table(glist), true, ...
            'Tmarkerlist','MarkerListTable');
    end

    function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx] = ismember(txt, glist);
        focalg = glist(idx);
    end

    function in_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end

    % function startDrag(hFig, evnt, targetAx)
    %     if isequal(evnt.Axes, targetAx)
    %         disp('Rotation started on target axes');
    %         hFig.WindowButtonMotionFcn = @(src,~) duringDrag(targetAx);
    %     end
    % end
    
    % function stopDragold(hFig, evnt, targetAx)
    %     %if isequal(evnt.Axes, targetAx)
    %         disp('Rotation stopped on target axes');
    %     %    hFig.WindowButtonMotionFcn = '';
    %     %end
    %     camPos = targetAx.CameraPosition;
    %     fprintf('Camera (%.2f, %.2f, %.2f)\n', camPos);
    %     [aa, bb]=view(targetAx);
    % end

%  copyobj(allchild(axOld(j)), axNew);


    function stopDrag(evnt, trackedAxes)
        thisAx = evnt.Axes;      % the axes that rotate3d thinks is active
        idxa = find(cellfun(@(h) isequal(h,thisAx), num2cell(trackedAxes)));       

        if ~isempty(idxa)
            if idxa > n
                tag = 'right';
                idxa = idxa - n;
            else
                tag = 'left';                
            end
            
            % fprintf('Rotation stopped on subplot %s #%d\n', tag, idxa);

            % trackedAxes(idxa).Tag
            % camPos = trackedAxes(idxa).CameraPosition;
            % [az,el] = view(trackedAxes(idxa));
            % fprintf('Camera (%.2f, %.2f, %.2f), View(%.2f,%.2f)\n', camPos, az, el);
            [yy, id]=ismember(tag, {'left','right'});
            assert(yy);
            [az, el] = view(ax{idxa, id});
            figure(hFig);
            if ~strcmp('Yes', gui.myQuestdlg(hFig, sprintf("Apply the same rotation to all tabs (%s plot)?", ...
                    tag))), return; end
               for kx = 1:n
                   if kx == idxa, continue; end
                   view(ax{kx, id}, [az, el]);
               end
        else
            disp('Rotation stopped on unknown axes (not in tracked list)');
        end
    end



    % function stopDrag_x(hFig, evnt)
    % 
    %     disp(evnt.Axes);           % see which axes handle MATLAB is giving you
    %     disp(evnt.Axes.Parent);    % see the parent (tab?)
    %     disp(evnt.Axes.Tag);       % check if tag survived        
    % 
    %     % axx = evnt.Axes;
    %     % if isempty(axx.Tag)
    %     %     disp('Rotation stopped on an axes without a Tag');
    %     % else
    %     %     disp(['Rotation stopped on ' axx.Tag]);
    %     % end
    % 
    %     % if isfield(evnt,'Axes') && strcmp(evnt.Axes.Tag, targetTag)
    %     % 
    %     %     targetTag
    %     % 
    %     %     disp(['Rotation stopped on 2 ' targetTag]);
    %     % 
    %     %     %camPos = evnt.Axes.CameraPosition;
    %     %     %fprintf('Camera (%.2f, %.2f, %.2f)\n', camPos);
    %     % 
    %     %     [az,el] = view(evnt.Axes);
    %     %     %fprintf('View az = %.2f, el = %.2f\n', az, el);
    %     %     % cleanup
    %     %     hFig.WindowButtonMotionFcn = '';
    %     %     figure(hFig);
    %     %     if strcmp('Yes', gui.myQuestdlg(hFig, "Apply to all tabs?"))
    %     %        for kx = 1:n
    %     %            view(ax{kx, 2}, [az, el]);
    %     %        end
    %     %     end
    %     % end
    % end
    
    % function duringDrag(ax)
    %     % Runs only while dragging on the chosen axes
    %     camPos = ax.CameraPosition;
    %     fprintf('Camera (%.2f, %.2f, %.2f)\n', camPos);
    % end

end




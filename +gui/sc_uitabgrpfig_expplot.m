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

pw1 = fileparts(mfilename('fullpath'));
%pth = fullfile(pw1, '..', 'resources', 'Misc', 'myTemplate.pptx');


hx=gui.myFigure;
hFig=hx.FigureHandle;
hFig.Position(3) = hFig.Position(3) * 1.8;

n = length(glist);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();
tab = cell(n,1);
ax0 = cell(n,1);
ax = cell(n,2);

idx = 1;
focalg = glist(idx);

% y = cell(n,1);
% for k=1:n
%     y{k} = sce.X(sce.g == glist(k), :);
% end

for k=1:n
    c = y{k};
    if issparse(c), c = full(c); end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',glist(k)));
    
    %{
    t = tiledlayout(1,2,'Parent',tab{k});
    ax1 = nexttile;
    hpl{k,1} = scatter3(s(:,1), s(:,2), s(:,3), 5, c, 'filled','Parent', ax1);
    ax2 = nexttile;
    hpl{k,2} = scatter(s(:,1), s(:,2), 5, c, 'filled','Parent', ax2);
    %}
    
    ax0{k} = axes('parent',tab{k});
    ax{k,1} = subplot(1,2,1);
    if size(s,2)>2
        scatter3(s(:,1), s(:,2), s(:,3), 5, c, 'filled');
    else
        scatter(s(:,1), s(:,2), 5, c, 'filled');
    end
    if ~isempty(cazcel)
        view(ax{k,1},cazcel(1),cazcel(2));
    end

    ax{k,2} = subplot(1,2,2);
        scatter(s(:,1), s(:,2), 5, c, 'filled');
        stem3(s(:,1), s(:,2), c, 'marker', 'none', 'color', 'm');
        hold on;
        scatter3(s(:,1), s(:,2), zeros(size(s(:,2))), 5, c, 'filled');
        
        title(ax{k,1}, glist(k));
        subtitle(ax{k,1}, gui.i_getsubtitle(c));
        title(ax{k,2}, glist(k));
        subtitle(ax{k,2}, gui.i_getsubtitle(c));

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


    gui.i_setautumncolor(c, a, true, any(c==0));
end
  
tabgp.SelectionChangedFcn=@displaySelection;

% hx.addCustomButton('off', [], "IMG00107.GIF", " ");
% hx.addCustomButton('off', @i_linksubplots, 'keyframes-minus.jpg', 'Link subplots');
hx.addCustomButton('off',  @in_genecards, 'www.jpg', 'GeneCards...');
%hx.addCustomButton('on', {@in_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
% hx.addCustomButton('off', @i_RescaleExpr, 'IMG00074.GIF', 'Rescale expression level [log2(x+1)]');
% hx.addCustomButton('off', @i_ResetExpr, 'plotpicker-geobubble2.gif', 'Reset expression level');
% hx.addCustomButton('off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');

hx.addCustomButton('off', @in_savedata, "floppy-disk-arrow-in.jpg", 'Save Gene List...');
%hx.addCustomButton('off', @in_savemainfig, "powerpoint.gif", 'Save Figure to PowerPoint File...');
%hx.addCustomButton('off', @in_savemainfigx, "File-Jpg--Streamline-Core-Gradient.png", 'Save Figure as Graphic File...');
%hx.addCustomButton('off', @in_set3dview, "tool_ellipse.gif", 'Set 3D View...');
hx.show(parentfig)

    % function in_set3dview(~, ~)
    %     [aa, bb] = view(ax{idx,2});
    %     answer = gui.myQuestdlg(hFig, 'Apply current view (azimuth and elevation angles) to all tabs?','');
    %     if strcmp(answer, 'Yes')
    %         for kx = 1:length(glist)
    %             if kx ~= idx
    %                 view(ax{kx,2}, [aa, bb]);
    %             end
    %         end
    %     end
    % end

    function in_savedata(~,~)
        gui.i_exporttable(table(glist), true, ...
            'Tmarkerlist','MarkerListTable');    
    end

    % function in_savemainfigx(~,~)
    %     answer = gui.myQuestdlg(hFig, 'Select Sub-plot to export:','', ...
    %         {'Left','Right','Cancel'},'Left');
    %     switch answer
    %         case 'Left'
    %             p = 1;
    %         case 'Right'
    %             p = 2;
    %         otherwise
    %             return;
    %     end
    %     [~,idx]=ismember(focalg, glist);     
    %     filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
    %     [filename, filepath] = uiputfile(filter,'Save Feature Plot', ...
    %         sprintf('FeaturePlot_%s', focalg));
    %     if ischar(filename)
    %         exportgraphics(ax{idx,p}, [filepath, filename]);
    %     end
    % end

    % function in_savemainfig(~,~)
    %     answer = gui.myQuestdlg(hFig, 'Export to PowerPoint?');
    %     if ~strcmp(answer,'Yes'), return; end
    % 
    %     fw=gui.gui_waitbar_adv;
    %         OUTppt = [tempname, '.pptx'];
    %         ppt = Presentation(OUTppt, pth);
    %         open(ppt);
    %         images=cell(n,1);
    %         warning off
    %     for kx=1:n
    %         gui.gui_waitbar_adv(fw,kx./n,"Processing "+glist(kx)+" ...");
    %         images{kx} = [tempname, '.png'];
    %         tabgp.SelectedTab=tab{kx};
    %         saveas(tab{kx},images{kx});
    %         slide3 = add(ppt, 'Small Title and Content');
    %         replace(slide3, 'Title', glist(kx));
    %         replace(slide3, 'Content', Picture(images{kx}));        
    %     end
    %         close(ppt);
    %         rptview(ppt);      
    %         gui.gui_waitbar_adv(fw);
    % end

    % function i_linksubplots(~,~)        
    %     hlink = linkprop([ax{idx,1},ax{idx,2}],{'CameraPosition','CameraUpVector'});
    % end

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

end


% function in_PickColorMap(~, ~, c)
%         list = {'parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', ...
%             'summer', 'autumn (default)', ...
%             'winter', 'jet'};
%         [indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
%             'PromptString', 'Select a colormap:', 'ListSize', [220, 300]);
%         if tf == 1
%             a = list{indx};
%             if strcmp(a, 'autumn (default)')
%                 a = 'autumn';
%             end
%             gui.i_setautumncolor(c, a);
%             setpref('scgeatoolbox', 'prefcolormapname', a);
%         end
%     end
% 

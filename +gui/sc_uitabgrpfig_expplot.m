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
%pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');


hx=gui.myFigure;
hFig=hx.FigHandle;
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
        view(ax{k,1}, [cazcel(1), cazcel(2)]);
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

end




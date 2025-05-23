function [hFig] = sc_multigroupings(sce, cx1, cx2, ttl1, ttl2, parentfig)

if nargin < 6, parentfig = []; end
if nargin < 5, ttl2 = ""; end
if nargin < 4, ttl1 = ""; end

c1 = cx1.c;
cL1 = cx1.cL;

c2 = cx2.c;
cL2 = cx2.cL;

hx=gui.myFigure;
hFig=hx.FigHandle;

ax1 = subplot(1, 2, 1);
h1 = gui.i_gscatter3(sce.s, c1, 1, 1);
if ~isempty(ttl1), title(ax1, strrep(ttl1,'_','\_')); end

%dt = datacursormode(hFig);
%dt.UpdateFcn = {@i_myupdatefcnx12};

ax2 = subplot(1, 2, 2);
h2 = gui.i_gscatter3(sce.s, c2, 1, 1);
if ~isempty(ttl2), title(ax2, strrep(ttl2,'_','\_')); end

dt = datacursormode(hFig);
dt.UpdateFcn = {@i_myupdatefcnx12};
sgtitle(sce.title);


kc1 = numel(unique(c1));
colormap(ax1, pkg.i_mycolorlines(kc1));
kc2 = numel(unique(c2));
colormap(ax2, pkg.i_mycolorlines(kc2));

evalin('base', 'h=findobj(gcf,''type'',''axes'');');
evalin('base', 'hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');

% h=findobj(hFig,'type','axes');
% linkprop(h,{'CameraPosition','CameraUpVector'});

rotate3d(hFig, 'on');
hFig.Position(3) = hFig.Position(3) * 2;
hBr = brush(hFig);
hBr.ActionPostCallback = {@onBrushAction, h1, h2};
hx.addCustomButton('off', @i_showclustlabel, "plotpicker-scatter.gif", "Show cluster labels");

hx.show(parentfig);


    function [txt] = i_myupdatefcnx12(Target, event_obj)
        % pos = event_obj.Position;
        if Target.Parent == ax1
            idx = event_obj.DataIndex;
            txt = cL1(c1(idx));
        elseif Target.Parent == ax2
            idx = event_obj.DataIndex;
            txt = cL2(c2(idx));
        end
    end

        function i_showclustlabel(~, ~)
            dtp1 = findobj(h1, 'Type', 'datatip');
            dtp2 = findobj(h2, 'Type', 'datatip');
            if ~isempty(dtp1) || ~isempty(dtp2)
                delete(dtp1);
                delete(dtp2);
                return;
            end
            if max(c1) < 50
                h1.DataTipTemplate.DataTipRows = dataTipTextRow('', cL1(c1));
                for i = 1:max(c1)
                    idx = find(c1 == i);
                    siv = sce.s(idx, :);
                    si = mean(siv, 1);
                    [k] = dsearchn(siv, si);
                    datatip(h1, 'DataIndex', idx(k));
                end
            end

            if max(c2) < 50
                h2.DataTipTemplate.DataTipRows = dataTipTextRow('', cL2(c2));
                for i = 1:max(c2)
                    idx = find(c2 == i);
                    siv = sce.s(idx, :);
                    si = mean(siv, 1);
                    [k] = dsearchn(siv, si);
                    datatip(h2, 'DataIndex', idx(k));
                end
            end

    end


    end

        function onBrushAction(~, event, h1, h2)
        if isequal(event.Axes.Children, h1)
            h2.BrushData = h1.BrushData;
        elseif isequal(event.Axes.Children, h2)
            h1.BrushData = h2.BrushData;
        end
        % if isprop(h1,'BrushData') && any(h1.BrushData)
        %     ptsSelected=logical(h1.BrushData);
        %     sum(ptsSelected)
        % end
    end

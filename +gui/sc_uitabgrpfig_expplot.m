function sc_uitabgrpfig_expplot(y, glist, s, parentfig, cazcel)

if nargin < 5, cazcel = []; end
if nargin < 4, parentfig = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

% https://www.mathworks.com/help/rptgen/ug/compile-a-presentation-program.html
if (ismcc || isdeployed) && pkg.i_isreportgenavailable('ppt'), makePPTCompilable(); end


hx = gui.myFigure(parentfig, true);

hFig = hx.FigHandle;
hFig.Position(3) = hFig.Position(3) * 1.8;

n = length(glist);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();
tab = cell(n,1);
ax = cell(n,2);

idx = 1;
focalg = glist(idx);

for k = 1:n
    c = y{k};
    if issparse(c), c = full(c); end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s', glist(k)));

    ax{k, 1} = subplot(1, 2, 1,'Parent', tab{k});

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

end
dorotation = false;
% Tabs added by "Show on the same figure...", one per side of the tabs.
mergedtab = gobjects(0);

trackedAxes = [ax{:,1}, ax{:,2}];   % only the right ones, or ax(:) if you want all
hRotate.ActionPostCallback = @(src,evnt) stopDrag(evnt, trackedAxes);

tabgp.SelectionChangedFcn=@displaySelection;
hx.addCustomButton('off', @in_genecards, 'www.jpg', 'GeneCards...');
hx.addCustomButton('off', @in_proteinstructure, 'hexagon_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Protein Structure...');
hx.addCustomButton('off', @in_savedata, "floppy-disk-arrow-in.jpg", 'Save Gene List...');
hx.addCustomButton('off', {@gui.callback_RunGeneAgent, glist}, "mw-microprocessor.jpg", 'Run GeneAgent...');
hx.addCustomButton('off', @in_mergetabs, 'Brightness-3--Streamline-Core.jpg', 'Show on the same figure...');
hx.addCustomButton('off', @in_colormap, 'www.jpg', 'Colormap...');
hx.addCustomButton('off', @in_syncrotation, 'www.jpg', 'Synchronize rotation angle...');
hx.show(parentfig)

function in_syncrotation(~, ~)
        dorotation = ~dorotation;
    end

function in_colormap(~, ~)
        gui.callback_PickColorMap(hFig, length(unique(c)), true, true);
    end

function in_mergetabs(~, ~)
        % New tabs in this window, not new windows. Rebuilt on every click,
        % so they show the gene tabs as they are now.
        delete(mergedtab(isvalid(mergedtab)));
        sides = ["All genes", "All genes, stem"];
        mergedtab = gobjects(1, 2);
        for kside = 1:2
            mergedtab(kside) = uitab(tabgp, 'Title', sides(kside));
            tl = tiledlayout(mergedtab(kside), 'flow');
            for kx = 1:n
                gui.i_cloneaxes(ax{kx, kside}, nexttile(tl));
            end
        end
        tabgp.SelectedTab = mergedtab(1);
    end

function in_savedata(~,~)
        gui.i_exporttable(table(glist), true, ...
            'Tmarkerlist', 'MarkerListTable', [], [], hFig);
    end

function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        [~,idx] = ismember(txt, glist);
        % A merged tab is not one gene; keep the last gene selected.
        if idx > 0, focalg = glist(idx); end
    end

function in_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end

function in_proteinstructure(~, ~)
        gui.i_viewprotein(focalg, ParentFig=hFig);
    end


function stopDrag(evnt, trackedAxes)
        if n < 2, return; end
        if ~dorotation, return; end
        thisAx = evnt.Axes;      % the axes that rotate3d thinks is active
        idxa = find(cellfun(@(h) isequal(h,thisAx), num2cell(trackedAxes)));

        if ~isempty(idxa)
            if idxa > n
                tag = 'right';
                idxa = idxa - n;
            else
                tag = 'left';
            end

            % trackedAxes(idxa).Tag
            % camPos = trackedAxes(idxa).CameraPosition;
            % [az,el] = view(trackedAxes(idxa));
            % fprintf('Camera (%.2f, %.2f, %.2f), View(%.2f,%.2f)\n', camPos, az, el);
            [yy, id]=ismember(tag, {'left','right'});
            assert(yy);
            [az, el] = view(ax{idxa, id});
            figure(hFig);
            if ~strcmp('Yes', gui.myQuestdlg(hFig, sprintf("Apply the same rotation to all tabs (% s plot)?", ...
                    tag))), return; end
               for kx = 1:n
                   if kx == idxa, continue; end
                   view(ax{kx, id}, [az, el]);
               end
        else
            disp('Rotation stopped on unknown axes (not in tracked list)');
        end
    end

end

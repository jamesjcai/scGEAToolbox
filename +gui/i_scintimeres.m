function i_scintimeres(s, g)

assert(size(s, 1) == length(g))
n = length(g);
FigureHandle = figure();
hAx = axes('Parent', FigureHandle);
c = ones(size(s, 1), 1);
[h] = gui.i_gscatter3(s, c, 1, 1, hAx);

kc = numel(unique(c));
colormap(pkg.i_mycolorlines(kc));

highlightedc = [];

dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

tb = uitoolbar('Parent', FigureHandle);
pkg.i_addbutton2fig(tb, 0, @refreshc, 'aa.gif', 'Clustering');
pkg.i_addbutton2fig(tb, 0, @highlightmax, 'aa.gif', 'Highlight top clusters');
pkg.i_addbutton2fig(tb, 0, @saveout, 'export.gif', 'Save top clusters');
pkg.i_addbutton2fig(tb, 0, @saveoutlist, 'export.gif', 'Save sorted gene list');

    function saveoutlist(~, ~)
        [~, idx] = sort(sum(s.^2, 2), 'descend');
        tmpName = [tempname, '.txt'];
        fid = fopen(tmpName, 'w');
        fprintf(fid, '%s\n', g(idx));
        fclose(fid);
        [status] = system(['notepad "', tmpName, '" &']);
        if status ~= 0
            if ~(ismcc || isdeployed)
                edit(tmpName);
            end
        end
end

        function saveout(~, ~)
            if isempty(highlightedc)
                errordlg('No selected cluters.');
                return;
            end
            tmpName = [tempname, '.txt'];
            fid = fopen(tmpName, 'w');
            for k = 1:length(highlightedc)
                a = sprintf('%s ', g(c == highlightedc(k)));
                fprintf(fid, 'Cluster %d\t%s\n', k, a);
            end
            %fprintf(fid,'%s',txt);
            fclose(fid);
            [status] = system(['notepad "', tmpName, '" &']);
            if status ~= 0
                edit(tmpName);
            end
    end


            function refreshc(~, ~)
                kcn = round(n/30);
                kc = gui.i_inputnumk(kcn, 10, 1000);
                if isempty(kc), return; end
                c = sc_cluster_s(s, kc, 'type', 'kmeans');
                h = gui.i_gscatter3(s, c, 1, hAx);
        end

                function highlightmax(~, ~)
                    kn = gui.i_inputnumk(10, 1, max(c));
                    a = grpstats(s, c, {@(x) mean(x, 1)});
                    d = sum(a.^2, 2);
                    [~, highlightedc] = sort(d, 'descend');
                    highlightedc = highlightedc(1:kn);
                    c2 = ones(size(c));
                    for k = 1:kn
                        c2(c == highlightedc(k)) = 2;
                    end
                    h = gui.i_gscatter3(s, c2, 1, hAx);
            end

                    function [txt] = i_myupdatefcnx(~, event_obj)
                        % pos = event_obj.Position;
                        idx = event_obj.DataIndex;
                        txt = g(idx);
                end

                end

function i_singlegraph(G1, figname, parentfig)
if nargin < 3, parentfig = []; end
if nargin < 2, figname = ''; end
if nargin < 1
    G1 = WattsStrogatz(100, 5, 0.15);
    G1.Nodes.Name = string((1:100)');
    G1.Edges.Weight = rand(size(G1.Edges, 1), 1) * 2;
end
import gui.*
import ten.*

%%

mfolder = fileparts(mfilename('fullpath'));

load(fullfile(mfolder, ...
    '../resources', 'tfome_tfgenes.mat'), 'tfgenes');


w = 3;
l = 1;

hFig = figure('Visible', 'off');
h1 = axes(hFig);
[p1] = drawnetwork(G1, h1);

%tb = uitoolbar(hFig);
tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');
pkg.i_addbutton2fig(tb, 'off', @ChangeFontSize, 'noun_font_size_591141.gif', 'ChangeFontSize');
pkg.i_addbutton2fig(tb, 'off', @ChangeWeight, 'noun_Weight_2243621.gif', 'ChangeWeight');
pkg.i_addbutton2fig(tb, 'off', @ChangeLayout, 'noun_Layout_792775.gif', 'ChangeLayout');
pkg.i_addbutton2fig(tb, 'off', @ChangeDirected, 'noun_directional_arrows_3497928.gif', 'ChangeDirected');
pkg.i_addbutton2fig(tb, 'off', @ChangeCutoff, 'noun_Pruners_2469297.gif', 'ChangeCutoff');
pkg.i_addbutton2fig(tb, 'off', @ChangeBox, 'noun_trim_3665385a.gif', 'Box off');
pkg.i_addbutton2fig(tb, 'off', @AnimateCutoff, 'noun_trim_3665385.gif', 'AnimateCutoff');
pkg.i_addbutton2fig(tb, 'off', @SaveAdj, 'export.gif', 'Export & save data');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');

title(h1,figname);

gui.i_movegui2parent(hFig, parentfig);

drawnow;
hFig.Visible=true;





    function SaveAdj(~, ~)
        if ~(ismcc || isdeployed)
            answer = questdlg('Export & save network to:', '', ...
                'Workspace', 'File', 'Workspace');
        else
            a = questdlg('Export & save network to file?');
            if strcmp(a, 'Yes')
                answer = 'File';
            else
                return;
            end
        end

        switch answer
            case 'Workspace'
                labels = {'Save adjacency matrix A1 to variable named:', ...
                    'Save graph G1 to variable named:', ...
                    'Save genelist g1 to variable named:'};
                A1 = full(adjacency(G1, 'weighted'));
                g1 = string(G1.Nodes.Name);
                vars = {'A', 'G', 'g'}; ...
                    values = {A1, G1, g1};
                msgfig = export2wsdlg(labels, vars, values);
                uiwait(msgfig);
            case 'File'
                [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
                if isequal(file, 0) || isequal(path, 0)
                    return;
                else
                    filename = fullfile(path, file);
                    pkg.m2cyto(G1, filename);
                end
            otherwise
                return;
        end
end


        function ChangeFontSize(~, ~)
            i_changefontsize(p1);
            %i_changefontsize(p2);
                function i_changefontsize(p)
                    if p.NodeFontSize >= 20
                        p.NodeFontSize = 7;
                    else
                        p.NodeFontSize = p.NodeFontSize + 1;
                    end
            end
        end

                function ChangeBox(~, ~)
                    if h1.Box
                        box(h1, 'off');
                        axis(h1, 'off');
                    else
                        box(h1, 'on');
                        axis(h1, 'on');
                    end
            end

                    function ChangeWeight(~, ~)
                        %a=3:10;
                        %w=a(randi(length(a),1));
                        w = w + 1;
                        if w > 10, w = 2; end
                        i_changeweight(p1, w);
                        %i_changeweight(p2,G2,w);
                            function i_changeweight(p, b)
                                %G.Edges.LWidths = abs(b*G.Edges.Weight/max(abs(G.Edges.Weight)));
                                %p.LineWidth = G.Edges.LWidths;
                                p.LineWidth = abs(b*p.LineWidth/max(abs(p.LineWidth)));
                        end
                    end

                            function ChangeLayout(~, ~)
                                a = ["auto", "layered", "subspace", "force", "circle"];
                                l = l + 1;
                                if l > 5, l = 1; end
                                switch a(l)
                                    case "force"
                                        p1.layout(a(l), 'Iterations', 250);
                                    otherwise
                                        p1.layout(a(l));
                                end
                                %p2.layout(a(l));
                                %p2.XData=p1.XData;
                                %p2.YData=p1.YData;
                                %p1.XData=p2.XData;
                                %p1.YData=p2.YData;
                        end

                                function ChangeDirected(~, ~)
                                    [p1, G1] = i_changedirected(p1, G1, h1);
                                    %[p2,G2]=i_changedirected(p2,G2,h2);
                                        function [p, G] = i_changedirected(p, G, h)
                                            x = p.XData;
                                            y = p.YData;
                                            if isa(G, 'digraph')
                                                A = adjacency(G, 'weighted');
                                                G = graph(0.5*(A + A.'), G.Nodes.Name);
                                                % p=plot(h,G);
                                                [p] = drawnetwork(G, h);
                                            end
                                            p.XData = x;
                                            p.YData = y;
                                    end
                                end

                                        function ChangeCutoff(~, ~)
                                            list = {'0.00 (show all edges)', ...
                                                '0.30', '0.35', '0.40', '0.45', ...
                                                '0.50', '0.55', '0.60', ...
                                                '0.65', '0.70', '0.75', '0.80', '0.85', ...
                                                '0.90', '0.95 (show 5% of edges)'};
                                            [indx, tf] = listdlg('ListString', list, ...
                                                'SelectionMode', 'single', ...
                                                'ListSize', [220, 300]);
                                            if tf
                                                if indx == 1
                                                    cutoff = 0;
                                                elseif indx == length(list)
                                                    cutoff = 0.95;
                                                else
                                                    cutoff = str2double(list(indx));
                                                end
                                                [p1] = i_replotg(p1, G1, h1, cutoff);
                                                %[p2]=i_replotg(p2,G2,h2,cutoff);
                                            end

                                    end


                                            function [p] = drawnetwork(G, h)
                                                p = plot(h, G);
                                                box off
                                                %layout(p,'force');
                                                %         if isa(G,'digraph')
                                                %             G.Nodes.NodeColors = outdegree(G)-indegree(G);
                                                %         else
                                                %             G.Nodes.NodeColors = degree(G);
                                                %         end
                                                %         p.NodeCData = G.Nodes.NodeColors;

                                                cc = repmat([0, 0.4470, 0.7410], G.numedges, 1);
                                                cc(G.Edges.Weight < 0, :) = repmat([0.8500, 0.3250, 0.0980], ...
                                                    sum(G.Edges.Weight < 0), 1);
                                                p.EdgeColor = cc;
                                                %       p.EdgeCData=ones(G.numedges,1);
                                                %       p.EdgeCData(G.Edges.Weight<0)=2;

                                                ix = ismember(string(upper(G.Nodes.Name)), tfgenes);
                                                if any(ix)
                                                    cc = repmat([0, 0, 0], G.numnodes, 1);
                                                    cc(ix, :) = repmat([1, 0, 0], sum(ix), 1);
                                                    p.NodeLabelColor = cc;
                                                end

                                                p.NodeFontSize = 2 * p.NodeFontSize;

                                                %title(h,sprintf('%d nodes',G.numnodes));
                                                % https://www.mathworks.com/matlabcentral/answers/296070-change-label-font-in-graph-plots
                                                %{
                                                nl = p.NodeLabel;
                                                p.NodeLabel = '';
                                                xd = get(p, 'XData');
                                                yd = get(p, 'YData');
                                                text(xd, yd, nl, 'FontSize',p.NodeFontSize,...
                                                    'FontWeight','bold',...
                                                    'HorizontalAlignment','left', ...
                                                    'VerticalAlignment','middle',...
                                                    'BackgroundColor','w','Margin',0.1);
                                                %}
                                                if ~isempty(G.Edges.Weight)
                                                    G.Edges.LWidths = abs(w*G.Edges.Weight/max(G.Edges.Weight));
                                                    p.LineWidth = G.Edges.LWidths;
                                                end
                                        end

                                                function AnimateCutoff(~, ~)
                                                    listc = 0.05:0.05:0.95;
                                                    f = waitbar(0, 'Cutoff = 0.05', 'Name', 'Edge Pruning...', ...
                                                        'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
                                                    setappdata(f, 'canceling', 0);

                                                    m = length(listc);
                                                    for k = 1:m
                                                        if getappdata(f, 'canceling')
                                                            break
                                                        end

                                                        cutoff = listc(k);
                                                        waitbar(k/m, f, sprintf('Cutoff = %g', cutoff));
                                                        try
                                                            p1 = i_replotg(p1, G1, h1, cutoff);
                                                        catch ME
                                                            disp(ME.message);
                                                        end
                                                        %p2=i_replotg(p2,G2,h2,cutoff);
                                                        pause(1);
                                                    end
                                                    %close(f)
                                                    delete(f)
                                            end

                                                    function [p, G] = i_replotg(p, G, h, cutoff)
                                                        a = h.Title.String;
                                                        x = p.XData;
                                                        y = p.YData;
                                                        A = adjacency(G, 'weighted');
                                                        A = ten.e_filtadjc(A, cutoff);
                                                        if issymmetric(A)
                                                            G = graph(A, G.Nodes.Name);
                                                        else
                                                            G = digraph(A, G.Nodes.Name);
                                                        end
                                                        % p=plot(h,G);
                                                        [p] = drawnetwork(G, h);
                                                        p.XData = x;
                                                        p.YData = y;
                                                        % h=gca;
                                                        title(a)
                                                end

                                                end


                                                    function h = WattsStrogatz(N, K, beta)
                                                    % H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
                                                    % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
                                                    %
                                                    % beta = 0 is a ring lattice, and beta = 1 is a random graph.

                                                    % Connect each node to its K next and previous neighbors. This constructs
                                                    % indices for a ring lattice.
                                                    s = repelem((1:N)', 1, K);
                                                    t = s + repmat(1:K, N, 1);
                                                    t = mod(t-1, N) + 1;

                                                    % Rewire the target node of each edge with probability beta
                                                    for source = 1:N
                                                        switchEdge = rand(K, 1) < beta;

                                                        newTargets = rand(N, 1);
                                                        newTargets(source) = 0;
                                                        newTargets(s(t == source)) = 0;
                                                        newTargets(t(source, ~switchEdge)) = 0;

                                                        [~, ind] = sort(newTargets, 'descend');
                                                        t(source, switchEdge) = ind(1:nnz(switchEdge));
                                                    end

                                                    h = graph(s, t);
                                                end

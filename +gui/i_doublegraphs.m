function [hFig] = i_doublegraphs(G1, G2, figname, parentfig)

if nargin < 4, parentfig = []; end
if nargin < 3, figname = ''; end
if nargin < 2
    G1 = WattsStrogatz(100, 5, 0.15);
    G2 = WattsStrogatz(100, 5, 0.15);
    G1.Nodes.Name = string((1:100)');
    G2.Nodes.Name = string((1:100)');
    G1.Edges.Weight = rand(size(G1.Edges, 1), 1) * 2;
    G2.Edges.Weight = rand(size(G2.Edges, 1), 1) * 2;
end
assert(isequal(G1.Nodes.Name, G2.Nodes.Name));
import gui.*
import ten.*

%%
mfolder = fileparts(mfilename('fullpath'));
load(fullfile(mfolder, ...
'..', 'assets', 'TFome', 'tfome_tfgenes.mat'), 'tfgenes');

w = 3;
l = 1;

hx=gui.myFigure(parentfig);
hFig = hx.FigHandle;
set(0, 'CurrentFigure', hFig);

tiledlayout(1, 2, 'TileSpacing', 'compact', ...
'Padding', 'compact')

h1 = nexttile;
[p1] = drawnetwork(G1, h1);

h2 = nexttile;
[p2] = drawnetwork(G2, h2);
p2.XData = p1.XData;
p2.YData = p1.YData;


hx.addCustomButton('off', @in_callback_ChangeFontSize, 'noun_font_size_591141.gif', 'Change Font Size of Nodes');
hx.addCustomButton('off', @in_callback_ChangeWeight, 'weight_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Width of Edges');
hx.addCustomButton('off', @in_callback_ChangeLayout, 'group_work_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Network Layout');
hx.addCustomButton('off', @in_callback_ChangeDirected, 'turn_sharp_right_17dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change to Undirected Network');
hx.addCustomButton('off', @in_callback_AnimateCutoff, 'movie.jpg', 'Animate Change and Select a Cutoff for Linked Edges');
hx.addCustomButton('off', @in_callback_ChangeCutoff, 'carpenter_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Cutoff to Trim Network');
hx.addCustomButton('off', @in_callback_SaveAdj, 'floppy-disk-arrow-in.jpg', 'Export & Save Data');
hx.addCustomButton('on',  @in_callback_RefreshAll, "refresh.jpg", "Refresh View");

hFig.Position(3) = hFig.Position(3) * 1.8;
if ~isempty(figname)
        figname = strrep(figname, '_', '\_');
        sgtitle(figname);
    end

hx.show(parentfig);

oldG1=[];
oldG2=[];
axistrig = true;


function in_callback_RefreshAll(~, ~)
        if ~isempty(oldG1)
            G1 = oldG1;
        end
        if ~isempty(oldG2)
            G2 = oldG2;
        end
        [p1] = drawnetwork(G1, h1);
        [p2] = drawnetwork(G2, h2);
        p2.XData = p1.XData;
        p2.YData = p1.YData;
    end


function in_callback_SaveAdj(~, ~)
        if ~(ismcc || isdeployed)
            labels = {'Save adjacency matrix A1 to variable named:', ...
                'Save adjacency matrix A2 to variable named:', ...
                'Save graph G1 to variable named:', ...
                'Save graph G2 to variable named:', ...
                'Save genelist g1 to variable named:', ...
                'Save genelist g2 to variable named:'};
            A1 = full(adjacency(G1, 'weighted'));
            A2 = full(adjacency(G2, 'weighted'));
            g1 = string(G1.Nodes.Name);
            g2 = string(G2.Nodes.Name);
            vars = {'A1', 'A2', 'G1', 'G2', 'g1', 'g2'}; ...
                values = {A1, A2, G1, G2, g1, g2};
            gui.i_export2wsdlg(hFig, labels, vars, values);
        else
            gui.myErrordlg(hFig, 'This function is not available for standalone application.');
        end
    end

function in_callback_ChangeFontSize(~, ~)
        i_changefontsize(p1);
        i_changefontsize(p2);
        function i_changefontsize(p)
            if p.NodeFontSize >= 20
                p.NodeFontSize = 7;
            else
                p.NodeFontSize = p.NodeFontSize + 1;
            end
        end
    end

function in_callback_ChangeWeight(~, ~)
        w = w + 1;
        if w > 10, w = 2; end
        i_changeweight(p1, w);
        i_changeweight(p2, w);
        function i_changeweight(p, b)
            p.LineWidth = abs(b*p.LineWidth/max(abs(p.LineWidth)));
        end
    end


function in_callback_ChangeLayout(~, ~)
        a = ["auto", "layered", "subspace", "force", "circle"];
        l = l + 1;
        if l > 5, l = 1; end
        switch a(l)
            case "force"
                p1.layout(a(l), 'Iterations', 500, ...
                    'WeightEffect', 'none', ...
                    'UseGravity', false);

            otherwise
                p1.layout(a(l));
        end

        p2.XData = p1.XData;
        p2.YData = p1.YData;
        p1.XData = p2.XData;
        p1.YData = p2.YData;
    end

function in_callback_ChangeDirected(~, ~)
        a1 = h1.Title.String;
        a2 = h2.Title.String;
        [p1, G1] = i_changedirected(p1, G1, h1);
        [p2, G2] = i_changedirected(p2, G2, h2);
        h1.Title.String = a1;
        h2.Title.String = a2;
        function [p, G] = i_changedirected(p, G, h)
            x = p.XData;
            y = p.YData;
            if isa(G, 'digraph')
                A = adjacency(G, 'weighted');
                G = graph(0.5*(A + A.'), G.Nodes.Name);
                [p] = drawnetwork(G, h);
            end
            p.XData = x;
            p.YData = y;
        end
    end

function in_callback_ChangeCutoff(~, ~)
        a1 = h1.Title.String;
        a2 = h2.Title.String;
        list = {'0.00 (show all edges)', ...
            '0.30', '0.35', '0.40', '0.45', ...
            '0.50', '0.55', '0.60', ...
            '0.65', '0.70', '0.75', '0.80', '0.85', ...
            '0.90', '0.95 (show 5% of edges)'};
        if gui.i_isuifig(parentfig)
            [indx, tf] = gui.myListdlg(hFig, list, '', [], false);
        else
            [indx, tf] = listdlg('ListString', list, ...
                'SelectionMode', 'single', 'ListSize', [220, 300]);
        end
        if tf
            if indx == 1
                cutoff = 0;
            elseif indx == length(list)
                cutoff = 0.95;
            else
                cutoff = str2double(list(indx));
            end
            [p1] = i_replotg(p1, G1, h1, cutoff);
            [p2] = i_replotg(p2, G2, h2, cutoff);
        end
        h1.Title.String = a1;
        h2.Title.String = a2;
    end

function [p] = drawnetwork(G, h)
        p = plot(h, G);
        n = size(G.Edges, 1);
        cc = repmat([0, 0.4470, 0.7410], n, 1);
        cc(G.Edges.Weight < 0, :) = repmat([0.8500, 0.3250, 0.0980], ...
            sum(G.Edges.Weight < 0), 1);
        p.EdgeColor = cc;

        i = ismember(string(upper(G.Nodes.Name)), tfgenes);
        if any(i)
            cc = repmat([0, 0, 0], G.numnodes, 1);
            cc(i, :) = repmat([1, 0, 0], sum(i), 1);
        end
        p.NodeFontSize = 2 * p.NodeFontSize;

        G.Edges.LWidths = abs(w*G.Edges.Weight/max(G.Edges.Weight));
        p.LineWidth = G.Edges.LWidths;

    end

function in_callback_AnimateCutoff(~, ~)
        listc = 0.05:0.05:0.95;
        % pkg.progressbar
        f = waitbar(0, 'Cutoff = 0.05', 'Name', 'Edge Pruning...', ...
            'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
        setappdata(f, 'canceling', 0);
        % Closed on every exit: closing the network window during a pause
        % made the next redraw throw and left this bar, Cancel button and
        % all, on screen.
        closeBar = onCleanup(@() delete(f(isvalid(f))));

        m = length(listc);
        for k = 1:m
            if ~isvalid(f) || getappdata(f, 'canceling') || ~isvalid(hFig)
                break   % Cancel keeps the current cutoff: that is the "select"
            end

            cutoff = listc(k);
            waitbar(k/m, f, sprintf('Cutoff = %g', cutoff));
            p1 = i_replotg(p1, G1, h1, cutoff);
            p2 = i_replotg(p2, G2, h2, cutoff);
            pause(2);
        end
    end

function [p, G] = i_replotg(p, G, h, cutoff)
        if ismember('Weight', G.Edges.Properties.VariableNames)
            if length(unique(G.Edges.Weight)) > 1
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
                [p] = drawnetwork(G, h);
                p.XData = x;
                p.YData = y;
                h.Title.String = a;
            end
        end
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

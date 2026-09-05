function st = i_cccstyle(G)
%I_CCCSTYLE Visual style for a causalCCC network, read off the graph itself.
%
%   st = i_cccstyle(G)
%
%   G is a network from SC_CAUSALCCCNET. This turns the columns that
%   function puts on it - Nodes.Color/Group/Side, Edges.Type/Sign/Weight -
%   into the arrays a plot needs, so the interactive view
%   (GUI.I_CAUSALCCCGRAPH) and the print view (GUI.I_NETWORKVISCCC) cannot
%   drift apart on what a colour means.
%
%   Fields:
%     IsBridge    numedges-by-1 logical, true for a known ligand-receptor
%                 link rather than an inferred dependency
%     IsAnchor    numnodes-by-1 logical, true for a ligand or receptor
%     NodeColor   numnodes-by-3 RGB, from Nodes.Color
%     NodeSize    numnodes-by-1 GraphPlot marker size; anchors larger,
%                 since they are what the view is built around
%     EdgeColor   numedges-by-3 RGB: gold bridge, orange negative,
%                 grey otherwise
%     EdgeWidth   numedges-by-1, inferred edges scaled by Weight. Fine as
%                 GraphPlot hairlines.
%     PrintWidth  numedges-by-1, bridge or not and nothing else. Drawn as
%                 real lines, the weight-scaled widths become a mat dense
%                 enough to bury the labels a print view is built around.
%     Legend      struct array of label/color/marker/lw rows, covering only
%                 what G actually contains
%
%   See also GUI.I_CAUSALCCCGRAPH, GUI.I_NETWORKVISCCC, SC_CAUSALCCCNET.

need = {'Color', 'Group', 'Side'};
missing = need(~ismember(need, G.Nodes.Properties.VariableNames));
if ~isempty(missing)
    error('gui:i_cccstyle:NotACausalcccGraph', ...
        ['G.Nodes is missing the %s column(s). This styling only applies ', ...
        'to a network built by SC_CAUSALCCCNET.'], strjoin(missing, ', '));
end
if ~all(ismember({'Type', 'Sign'}, G.Edges.Properties.VariableNames))
    error('gui:i_cccstyle:NotACausalcccGraph', ...
        ['G.Edges is missing the Type or Sign column. This styling only ', ...
        'applies to a network built by SC_CAUSALCCCNET.']);
end

st.IsBridge = G.Edges.Type == "bridge";
st.IsAnchor = ismember(G.Nodes.Group, ["ligand", "receptor"]);

% ---- nodes -------------------------------------------------------------
st.NodeColor = i_hex2rgb(G.Nodes.Color);
st.NodeSize = repmat(5, G.numnodes, 1);
st.NodeSize(st.IsAnchor) = 9;

% ---- edges -------------------------------------------------------------
st.EdgeColor = repmat([0.55, 0.55, 0.55], G.numedges, 1);
neg = G.Edges.Sign < 0;
st.EdgeColor(neg, :) = repmat([0.85, 0.33, 0.10], sum(neg), 1);
st.EdgeColor(st.IsBridge, :) = repmat([0.85, 0.65, 0.05], sum(st.IsBridge), 1);

inferred = ~st.IsBridge;
st.EdgeWidth = ones(G.numedges, 1);
if any(inferred) && max(G.Edges.Weight(inferred)) > 0
    st.EdgeWidth(inferred) = rescale(G.Edges.Weight(inferred), 1, 4);
end
st.EdgeWidth(st.IsBridge) = 3;

st.PrintWidth = ones(G.numedges, 1);
st.PrintWidth(st.IsBridge) = 3;

% ---- legend ------------------------------------------------------------
% Only the rows G actually uses: a network with no negative edges, because
% the Method was unsigned or it simply has none, should not advertise a
% colour it never draws.
labels = struct( ...
    'ligand',        'ligand (sender anchor)', ...
    'sender_gene',   'sender-side gene', ...
    'receptor',      'receptor (receiver anchor)', ...
    'receiver_gene', 'receiver-side gene', ...
    'metadata',      'metadata pivot');
st.Legend = struct('label', {}, 'color', {}, 'marker', {}, 'lw', {});
for grp = unique(G.Nodes.Group, 'stable')'
    key = char(grp);
    if isfield(labels, key), lab = labels.(key); else, lab = key; end
    st.Legend(end+1) = struct('label', lab, 'marker', 'o', 'lw', 0.5, ...
        'color', st.NodeColor(find(G.Nodes.Group == grp, 1), :));
end
if any(st.IsBridge)
    st.Legend(end+1) = struct('label', 'ligand-receptor link (known)', ...
        'color', [0.85, 0.65, 0.05], 'marker', '-', 'lw', 3);
end
if any(inferred & ~neg)
    st.Legend(end+1) = struct('label', 'inferred dependency', ...
        'color', [0.55, 0.55, 0.55], 'marker', '-', 'lw', 1);
end
if any(inferred & neg)
    st.Legend(end+1) = struct('label', 'inferred, negative', ...
        'color', [0.85, 0.33, 0.10], 'marker', '-', 'lw', 1);
end
end


function rgb = i_hex2rgb(hexstr)
hexstr = string(hexstr);
rgb = zeros(numel(hexstr), 3);
for k = 1:numel(hexstr)
    h = char(erase(hexstr(k), "#"));
    rgb(k, :) = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))] / 255;
end
end

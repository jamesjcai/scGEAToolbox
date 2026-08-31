function [G, info] = e_sce2causalcccview(sce, senders, receivers, parentfig, varargin)
%GUI.E_SCE2CAUSALCCCVIEW  Build and plot a local causalCCC-style network.
%
%   G = GUI.E_SCE2CAUSALCCCVIEW(sce, senders, receivers)
%   G = GUI.E_SCE2CAUSALCCCVIEW(sce, senders, receivers, parentfig)
%   G = GUI.E_SCE2CAUSALCCCVIEW(sce, senders, receivers, parentfig, Name=Value)
%   [G, info] = GUI.E_SCE2CAUSALCCCVIEW(...)
%
%   GUI-facing wrapper around SC_CAUSALCCCNET: shows a waitbar while the
%   network is built, an error dialog if it fails, and otherwise plots it
%   with GUI.I_CAUSALCCCGRAPH. Every Name=Value pair SC_CAUSALCCCNET takes
%   (GroupBy, LigandReceptor, SenderGenes, ReceiverGenes, GeneSelection,
%   MaxGenes, MaxPairs, MinFrac, Metadata, MIBins, MIPermutations, Method,
%   TopKPerNode, EdgeCutoff) is accepted here and passed straight through.
%
%   parentfig embeds the plot in an existing GUI figure, as elsewhere in
%   this toolbox (e.g. SC_GRNVIEW); pass [] or omit it to open a new one.
%
%   G and info are empty if the build fails (the error dialog already told
%   the user why); check ~isempty(G) before using them in a script.
%
%   EXAMPLE
%     load('example_data/new_example_sce.mat', 'sce');
%     G = gui.e_sce2causalcccview(sce, "Macrophages", "Beta cells");
%
%   See also: SC_CAUSALCCCNET, PKG.SC_SCE2CAUSALCCC, GUI.I_CAUSALCCCGRAPH

if nargin < 4, parentfig = []; end

fw = gui.myWaitbar(parentfig);
try
    [G, info] = sc_causalcccnet(sce, senders, receivers, varargin{:});
catch ME
    gui.myWaitbar(parentfig, fw, true);
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    G = [];
    info = [];
    return
end
gui.myWaitbar(parentfig, fw);

ttl = sprintf('%s -> %s', strjoin(senders, '+'), strjoin(receivers, '+'));
subttl = sprintf('%d senders, %d receivers; %d nodes, %d edges', ...
    info.nSenderCells, info.nReceiverCells, G.numnodes, G.numedges);
gui.i_causalcccgraph(G, [string(ttl); string(subttl)], parentfig);
end

function callback_BuildGRNAllGenes(src, cellIdx)
% CALLBACK_BUILDGRNALLGENES  Build a genome-wide GRN from all genes.
%
% The all-genes branch of Network > Build Gene Regulatory Network (GRN)...:
% gui.callback_BuildGeneNetwork has already picked the cells (logical
% CELLIDX, one per cell; all cells when omitted). Here: pick the method,
% then one waitbar while the network is built. gui.i_showgrn then saves it
% with its gene list to the working folder and draws it as a sketch
% (sc_grnsketch); the sketch's toolbar exports it to the workspace and
% shows the saved file in its folder.
%
% Methods
%   PCR           net.pcrnet - the signed scTenifoldNet network used by the
%                 scTenifold comparison and knockout tools.
%   Denoised PCR  net.pcrnet_denoised - PCR networks from 10 cell
%                 subsamples, merged by CP tensor decomposition (the full
%                 scTenifoldNet construction; needs Tensor Toolbox).
%   IDS           net.idsnet - unsigned, nonlinear co-dependence network;
%                 not an scTenifold input, so its file also stores METHOD.
%
% See also gui.callback_BuildGeneNetwork, net.pcrnet, net.pcrnet_denoised,
%   net.idsnet, gui.i_showgrn, sc_grnsketch, gui.callback_scTenifoldNetView.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if nargin < 2 || isempty(cellIdx), cellIdx = true(sce.NumCells, 1); end

% --- Method
methods = {'PCR (signed; scTenifoldNet)', ...
    'Denoised PCR (10 subsamples + tensor decomposition)', ...
    'IDS (unsigned, nonlinear)'};
prompt = sprintf(['How should the network over all %d genes be built?\n\n', ...
    'PCR is fast and gives the signed network the scTenifold comparison ', ...
    'and knockout tools use. Denoised PCR builds PCR networks on 10 ', ...
    'subsamples of 500 cells and merges them; slower, needs Tensor ', ...
    'Toolbox. IDS scores unsigned, nonlinear dependence. ', ...
    '[scTenifoldNet, PMID:33336197]'], numel(sce.g));
[indx, tf] = gui.myListdlg(FigureHandle, methods, 'Network Method', ...
    1, false, true, [420, 160], prompt);
if tf ~= 1 || isempty(indx), return; end

if indx == 2 && ~i_checktensortoolbox(FigureHandle)
    return;
end

% --- Build
X = sce.X(:, cellIdx);
g = sce.g;
fw = gui.myWaitbar(FigureHandle);
try
    switch indx
        case 1
            disp('>> A = net.pcrnet(log1p(sc_norm(X)));')
            X = log1p(sc_norm(X));
            A = net.pcrnet(X, 3, false, true, false, false, pkg.i_usegpu(X));
            method = '';
            tag = 'pcr';
        case 2
            disp('>> A = net.pcrnet_denoised(X, ''savegrn'', false);')
            A = net.pcrnet_denoised(X, 'savegrn', false);
            method = '';
            tag = 'denoisedpcr';
        case 3
            % net.idsnet scales each gene into [0, 1] itself; library-size
            % normalised log1p data is the input it was tuned on.
            disp('>> A = net.idsnet(log1p(sc_norm(X)));')
            A = net.idsnet(log1p(sc_norm(X)));
            method = 'ids';
            tag = 'ids';
        otherwise
            % myListdlg returns an index into METHODS; nothing else can occur
            error('gui:callback_BuildGRNAllGenes:BadMethod', 'Unknown method.');
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

try
    gui.i_showgrn(A, g, method, tag, FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end
end


function ok = i_checktensortoolbox(FigureHandle)
ok = true;
try
    ten.check_tensor_toolbox;
catch
    gui.i_installtensortoolbox(FigureHandle);
    try
        ten.check_tensor_toolbox;
    catch ME
        gui.myErrordlg(FigureHandle, ME.message);
        ok = false;
    end
end
end


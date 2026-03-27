function A = grnformer(X, tf_idx, varargin)
%GRNFORMER  GRNFormer: Graph Transformer for Gene Regulatory Network inference
%
%  A = net.grnformer(X, tf_idx)
%  A = net.grnformer(X, tf_idx, 'GroundTruth', GT, ...)
%
%  Inputs
%    X        [G x C] genes-by-cells expression matrix (raw or log-counts)
%    tf_idx   integer indices (or logical vector) of TF genes within X rows
%
%  Name-Value options
%    'GroundTruth'    [G x G] logical/sparse known regulatory adjacency
%    'Weights'        pretrained params struct (skips training)
%    'SubgraphSize'   TF-Walker neighborhood size      (default 100)
%    'CorrThreshold'  GCEN Pearson |r| cutoff           (default 0.1)
%    'EmbedDim'       Gene-Transcoder embedding dim     (default 64)
%    'LatentDim'      VAE latent dim                    (default 16)
%    'NumHeads'       attention heads per TransConv      (default 4)
%    'EncBlocks'      encoder TransConv blocks           (default 4)
%    'DecBlocks'      decoder TransConv blocks           (default 3)
%    'Epochs'         training epochs                   (default 100)
%    'BatchSize'      subgraph mini-batch size           (default 8)
%    'LearningRate'   Adam learning rate                (default 0.001)
%    'MaxCells'       max cells used in Gene-Transcoder  (default 200)
%    'L1Weight'       L1 regularisation weight           (default 1e-5)
%    'Verbose'        print training progress            (default true)
%
%  Output
%    A        [G x G] sparse matrix of TF-to-gene regulatory scores
%
%  Reference: Hegde A & Cheng J. GRNFormer: Accurate Gene Regulatory
%    Network Inference Using Graph Transformer.
%    Bioinformatics 2026, doi: 10.1093/bioinformatics/btag144
%
%  See also: sc_grn, net.pearsonnet, net.genie3

% ---- argument parsing -----------------------------------------------
p = inputParser;
addRequired(p,'X');
addRequired(p,'tf_idx');
addParameter(p,'GroundTruth',  []);
addParameter(p,'Weights',      []);
addParameter(p,'SubgraphSize', 100);
addParameter(p,'CorrThreshold',0.1);
addParameter(p,'EmbedDim',     64);
addParameter(p,'LatentDim',    16);
addParameter(p,'NumHeads',     4);
addParameter(p,'EncBlocks',    4);
addParameter(p,'DecBlocks',    3);
addParameter(p,'Epochs',       100);
addParameter(p,'BatchSize',    8);
addParameter(p,'LearningRate', 1e-3);
addParameter(p,'MaxCells',     200);
addParameter(p,'L1Weight',     1e-5);
addParameter(p,'Verbose',      true);
parse(p, X, tf_idx, varargin{:});
o = p.Results;

X = double(full(X));
[G, C] = size(X);

% Build TF list
if islogical(tf_idx)
    tf_list = find(tf_idx(:));
else
    tf_list = unique(tf_idx(:));
end
tf_mask = false(G,1);
tf_mask(tf_list) = true;
T = numel(tf_list);
if T == 0, error('grnformer:noTFs','tf_idx identified zero TF genes.'); end

if o.Verbose
    fprintf('[GRNFormer] %d genes x %d cells, %d TFs\n', G, C, T);
end

% ---- 1. arcsinh normalisation ----------------------------------------
Xn = asinh(X);

% ---- 2. Gene co-expression network -----------------------------------
if o.Verbose, fprintf('[GRNFormer] Building GCEN (thr=%.2f)...\n', o.CorrThreshold); end
GCEN = buildGCEN(Xn, o.CorrThreshold);          % G x G sparse Pearson

% ---- 3. Cell subsampling for Gene-Transcoder -------------------------
Cuse = min(C, o.MaxCells);
rng_state = rng;                                 % reproducible cell subset
rng(42);
cell_sel = randperm(C, Cuse);
rng(rng_state);
Xc = Xn(:, cell_sel);                           % G x Cuse

% ---- 4. TF-Walker training subgraphs ---------------------------------
N_sub = o.SubgraphSize;
if o.Verbose, fprintf('[GRNFormer] TF-Walker sampling (%d TFs)...\n', T); end
[sg_nodes, sg_edges] = tfWalkerTrain(GCEN, tf_list, N_sub);

% ---- 5. Node feature extraction per subgraph -------------------------
n_sg = numel(sg_nodes);
node_feats = cell(n_sg,1);
for k = 1:n_sg
    nidx = sg_nodes{k};
    Xs = cellwiseZscore(Xc(nidx,:));            % N_sub x Cuse, z-scored
    tff = double(tf_mask(nidx));                % N_sub x 1  TF identity flag
    node_feats{k} = [Xs, tff];                  % N_sub x (Cuse+1)
end

% ---- 6. Ground-truth labels per subgraph -----------------------------
GT = o.GroundTruth;
if isempty(GT)
    if o.Verbose
        fprintf('[GRNFormer] No GroundTruth provided -- using top-20%% co-expression as proxy labels.\n');
    end
    Rfull = abs(full(GCEN));
    Rfull(1:G+1:end) = 0;             % no self-loops
    vals = Rfull(Rfull > 0);
    if isempty(vals)
        thr80 = o.CorrThreshold;
    else
        thr80 = quantile(vals, 0.80); % top 20% of co-expressed pairs
    end
    GTmat = double(Rfull >= thr80);
else
    GTmat = double(full(GT));
end

sg_labels = cell(n_sg,1);
for k = 1:n_sg
    nidx = sg_nodes{k};
    L = GTmat(nidx, nidx);
    % keep only TF-initiated edges (TF in row)
    L = L .* double(tf_mask(nidx));
    L(1:N_sub+1:end) = 0;                       % no self-loops
    sg_labels{k} = L;
end

% Discard subgraphs with no positive labels
ok = cellfun(@(l) any(l(:) > 0), sg_labels);
sg_nodes  = sg_nodes(ok);
sg_edges  = sg_edges(ok);
sg_labels = sg_labels(ok);
node_feats= node_feats(ok);
n_sg = numel(sg_nodes);
if n_sg == 0
    error('grnformer:noLabels', ...
        'No subgraphs have positive labels. Lower CorrThreshold or provide GroundTruth.');
end
if o.Verbose, fprintf('[GRNFormer] %d valid training subgraphs.\n', n_sg); end

% ---- 7. Model parameters ---------------------------------------------
d_in  = Cuse + 1;
d_emb = o.EmbedDim;
d_lat = o.LatentDim;
nh    = o.NumHeads;

if ~isempty(o.Weights)
    params = o.Weights;
    do_train = false;
    if o.Verbose, fprintf('[GRNFormer] Using pretrained weights.\n'); end
else
    params = initParams(d_in, d_emb, d_lat, nh, o.EncBlocks, o.DecBlocks);
    do_train = true;
end

% ---- 8. Training -------------------------------------------------------
if do_train
    params = trainModel(params, node_feats, sg_edges, sg_labels, ...
        N_sub, d_emb, d_lat, nh, o.EncBlocks, o.DecBlocks, ...
        o.Epochs, o.BatchSize, o.LearningRate, o.L1Weight, o.Verbose);
end

% ---- 9. Inference -------------------------------------------------------
if o.Verbose, fprintf('[GRNFormer] Running inference...\n'); end
[inf_nodes, inf_edges] = tfWalkerInfer(GCEN, tf_list, N_sub);

score_sum = zeros(G,G,'single');
score_cnt = zeros(G,G,'uint16');

for k = 1:numel(inf_nodes)
    nidx = inf_nodes{k};
    E    = inf_edges{k};
    nf   = [cellwiseZscore(Xc(nidx,:)), double(tf_mask(nidx))];

    probs = forwardInfer(params, nf, E, d_emb, d_lat, nh, o.EncBlocks, o.DecBlocks);
    probs = max(0, min(1, probs));

    w = single(probs .* (probs > 0.3) .* double(tf_mask(nidx)));
    score_sum(nidx, nidx) = score_sum(nidx, nidx) + w;
    score_cnt(nidx, nidx) = score_cnt(nidx, nidx) + 1;
end

nz = score_cnt > 0;
S  = zeros(G,G,'double');
S(nz) = double(score_sum(nz)) ./ double(score_cnt(nz));
S = S .* double(tf_mask);                       % TF-initiated only
S(1:G+1:end) = 0;                               % no self-loops
mx = max(S(:)); mn = min(S(:));
if mx > mn, S = (S - mn) / (mx - mn); end

A = sparse(S);
if o.Verbose, fprintf('[GRNFormer] Done. Non-zero edges: %d\n', nnz(A)); end
end  % grnformer


% ======================================================================
%  GCEN  /  TF-Walker
% ======================================================================

function GCEN = buildGCEN(X, thr)
% Pearson correlation, |r| > thr, no self-loops → sparse G x G matrix
[G, C] = size(X);
mu  = mean(X,2);
sg  = std(X,0,2);  sg(sg<1e-10) = 1;
Xz  = (X - mu) ./ sg;
R   = (Xz * Xz') / max(C-1, 1);
R   = max(min(R,1),-1);
R(abs(R) < thr) = 0;
R(1:G+1:end) = 0;
GCEN = sparse(R);
end

% -----------------------------------------------------------------------
function [sg_nodes, sg_edges] = tfWalkerTrain(GCEN, tf_list, N_sub)
% Training version: random selection when neighbours > N_sub-1
T   = numel(tf_list);
adj = logical(abs(GCEN) > 0);
sg_nodes = cell(T,1);
sg_edges = cell(T,1);

for ti = 1:T
    tf  = tf_list(ti);
    nb  = expandNeighbors(adj, tf, N_sub-1);
    if numel(nb) > N_sub-1
        nb = nb(randperm(numel(nb), N_sub-1));
    end
    nidx = padSubgraph([tf, nb(:)'], N_sub);
    sg_nodes{ti} = nidx;
    sg_edges{ti} = full(GCEN(nidx, nidx));
end
end

% -----------------------------------------------------------------------
function [sg_nodes, sg_edges] = tfWalkerInfer(GCEN, tf_list, N_sub)
% Inference version: deterministic, sequential chunking of all neighbours
adj = logical(abs(GCEN) > 0);
sg_nodes = {};
sg_edges = {};

for ti = 1:numel(tf_list)
    tf = tf_list(ti);
    nb = find(adj(tf,:));
    if isempty(nb), continue; end

    all_nb = [tf, nb(:)'];
    % Slide a window of N_sub over the full neighbourhood list
    for s = 1:N_sub-1:numel(all_nb)
        chunk = all_nb(s : min(s+N_sub-2, numel(all_nb)));
        nidx  = padSubgraph(unique([tf, chunk],'stable'), N_sub);
        sg_nodes{end+1} = nidx; %#ok<AGROW>
        sg_edges{end+1} = full(GCEN(nidx, nidx)); %#ok<AGROW>
    end
end
end

% -----------------------------------------------------------------------
function nb = expandNeighbors(adj, seed, target)
% BFS from seed until at least 'target' neighbours found (max 6 hops)
nb = [];
frontier = seed;
visited  = false(size(adj,1),1);
visited(seed) = true;
for hop = 1:6
    new_nb = find(any(adj(frontier,:), 1));
    new_nb = new_nb(~visited(new_nb));
    if isempty(new_nb), break; end
    visited(new_nb) = true;
    nb = [nb, new_nb]; %#ok<AGROW>
    frontier = new_nb;
    if numel(nb) >= target, break; end
end
end

% -----------------------------------------------------------------------
function nidx = padSubgraph(nidx, N)
% Pad (by repeating last entry) or trim to exactly N nodes
if numel(nidx) < N
    nidx = [nidx, repmat(nidx(end), 1, N-numel(nidx))];
else
    nidx = nidx(1:N);
end
end

% -----------------------------------------------------------------------
function Xz = cellwiseZscore(X)
% Z-score across genes for each cell (column-wise standardisation)
mu = mean(X,1);
sg = std(X,0,1);  sg(sg<1e-10) = 1;
Xz = (X - mu) ./ sg;
end


% ======================================================================
%  PARAMETER INITIALISATION
% ======================================================================

function params = initParams(~, d_emb, d_lat, nh, n_enc, n_dec)
% He initialisation for all learnable weight matrices

he = @(m,n) dlarray(randn(m,n) * sqrt(2/max(m,1)));
zb = @(n)   dlarray(zeros(1,n));
on = @(n)   dlarray(ones(1,n));

% ---- Gene-Transcoder ---------------------------------------------------
% Conv1D kernel=1: each cell's scalar expression → d_emb channels
params.gt.W_conv  = he(1, d_emb);             % 1 input channel → d_emb
params.gt.b_conv  = zb(d_emb);
% single transformer encoder layer (self-attention over cells per gene)
params.gt.Wq = he(d_emb, d_emb);
params.gt.Wk = he(d_emb, d_emb);
params.gt.Wv = he(d_emb, d_emb);
params.gt.Wo = he(d_emb, d_emb);
params.gt.ffW1 = he(d_emb, d_emb*2);  params.gt.ffb1 = zb(d_emb*2);
params.gt.ffW2 = he(d_emb*2, d_emb);  params.gt.ffb2 = zb(d_emb);
params.gt.ln1g = on(d_emb);  params.gt.ln1b = zb(d_emb);
params.gt.ln2g = on(d_emb);  params.gt.ln2b = zb(d_emb);

% ---- GraViTAE encoder blocks ------------------------------------------
for bi = 1:n_enc
    params.enc(bi) = makeTransConvParams(d_emb, nh, he, zb, on);
end
params.enc_mu_n = he(d_emb, d_lat);    % node → latent mean
params.enc_lv_n = he(d_emb, d_lat);    % node → latent log-var
params.enc_mu_e = he(2*d_emb, d_lat);  % edge (concat node pair) → mean
params.enc_lv_e = he(2*d_emb, d_lat);  % edge → log-var

% ---- GraViTAE decoder blocks ------------------------------------------
for bi = 1:n_dec
    params.dec(bi) = makeTransConvParams(d_lat, nh, he, zb, on);
end
params.dec_W1 = he(d_lat, d_emb);  params.dec_b1 = zb(d_emb);
params.dec_W2 = he(d_emb, d_emb);  params.dec_b2 = zb(d_emb);
end

% -----------------------------------------------------------------------
function bp = makeTransConvParams(d, ~, he, zb, on)
% Create parameter block for one TransConv block of width d
bp.Wq    = he(d,d);  bp.Wk  = he(d,d);  bp.Wv = he(d,d);
bp.We_k  = he(1,d);                      % edge → key contribution
bp.We_v  = he(1,d);                      % edge → value contribution
bp.Wo    = he(d,d);
bp.ffW1  = he(d,d*2);  bp.ffb1 = zb(d*2);
bp.ffW2  = he(d*2,d);  bp.ffb2 = zb(d);
bp.ln1g  = on(d);  bp.ln1b = zb(d);
bp.ln2g  = on(d);  bp.ln2b = zb(d);
end


% ======================================================================
%  FORWARD PASS
% ======================================================================

function Z_emb = geneTranscoder(p, nf, d_emb, nh)
% nf:    N x (Cuse+1)  node feature matrix (expression + TF flag)
% Z_emb: N x d_emb     fixed gene embeddings
%
% Architecture: pointwise Conv1D → single Transformer layer → mean pool

% 1. Conv1D (kernel=1): treat each cell as a token, transform d_in → d_emb
%    nf is N x (Cuse+1); after conv each of the N genes has Cuse+1 tokens of d_emb
nf_dl = dlarray(nf);            % N x (Cuse+1)
% Equivalent: for each of the N genes, project each scalar to d_emb
% Vectorised: nf is N x (Cuse+1), treat as batch over cells
% Result: N x (Cuse+1) x d_emb via broadcast multiply
%   tokens[g, c, :] = nf[g,c] * W_conv + b_conv
% Reshape nf to (N*(Cuse+1)) x 1 then linear → (N*(Cuse+1)) x d_emb
nf_flat = reshape(nf_dl, [], 1);                    % N*(Cuse+1) x 1
tokens  = nf_flat * p.W_conv + p.b_conv;            % N*(Cuse+1) x d_emb
[N, T_cell] = size(nf);
tokens  = reshape(tokens, N, T_cell, d_emb);        % N x T_cell x d_emb

% 2. Single transformer layer: self-attention over T_cell tokens, per gene
% Process all N genes simultaneously by merging N into batch dimension
% Reshape: (N*T_cell) x d_emb, but attention is within each gene's T_cell tokens
% Use a block-diagonal / batched approach: loop over genes
% For N=100, T_cell<=200, this loop is N=100 iterations (manageable)

% Pre-allocate output as a dlarray-compatible accumulation
gene_embeds = zeros(N, d_emb, 'like', extractdata(tokens));
for g = 1:N
    tok_g = squeeze(tokens(g,:,:));   % T_cell x d_emb
    tok_g = dlarray(tok_g);           % re-wrap (needed inside dlfeval)
    % Self-attention
    Q = tok_g * p.Wq;                 % T x d
    K = tok_g * p.Wk;
    V = tok_g * p.Wv;
    dh = d_emb / nh;
    S = (Q * K') / sqrt(dh);         % T x T attention scores
    A = rowsoftmax(S);               % row-wise softmax
    H = A * V * p.Wo;                % T x d_emb

    % Residual + layer norm 1
    H1 = layernorm(tok_g + H, p.ln1g, p.ln1b);

    % FFN
    ff = leakyrelu(H1 * p.ffW1 + p.ffb1, 0.01);
    ff = ff * p.ffW2 + p.ffb2;

    % Residual + layer norm 2
    out_g = layernorm(H1 + ff, p.ln2g, p.ln2b);

    % Mean pool over cells → 1 x d_emb embedding for gene g
    gene_embeds(g,:) = extractdata(mean(out_g, 1));
end
Z_emb = dlarray(gene_embeds);   % N x d_emb
end

% -----------------------------------------------------------------------
function [mu_n, lv_n, mu_e, lv_e] = gravitaeEncode(p, Z, E, n_enc, nh)
% Z:    N x d_emb  node embeddings from Gene-Transcoder
% E:    N x N      edge features (Pearson correlations)
% Output: per-node and per-edge latent parameters

H = Z;
for bi = 1:n_enc
    H = transConvBlock(p.enc(bi), H, E, nh);
end

% Node latent
mu_n = H * p.enc_mu_n;              % N x d_lat
lv_n = H * p.enc_lv_n;              % N x d_lat

% Edge latent: per directed edge (i,j) = concat(H_i, H_j) → linear
N = size(H,1);

% Build pairwise concat: N x N x 2*d_emb (vectorised)
% Expand: H_i replicated across columns, H_j replicated across rows
Hi_exp = permute(repmat(H, 1, 1, N), [1,3,2]);   % N x N x d_emb
Hj_exp = permute(repmat(H, 1, 1, N), [3,1,2]);   % N x N x d_emb
Hpair  = cat(3, Hi_exp, Hj_exp);                  % N x N x 2*d_emb
Hp_flat = reshape(Hpair, N*N, 2*size(H,2));       % (N^2) x 2*d_emb

mu_e_flat = Hp_flat * p.enc_mu_e;                 % (N^2) x d_lat
lv_e_flat = Hp_flat * p.enc_lv_e;
d_lat = size(mu_e_flat,2);
mu_e = reshape(mu_e_flat, N, N, d_lat);           % N x N x d_lat
lv_e = reshape(lv_e_flat, N, N, d_lat);
end

% -----------------------------------------------------------------------
function [Z_prime, E_out] = gravitaeDecode(p, z_n, z_e_mean, E, n_dec, nh)
% z_n:      N x d_lat  sampled node latent
% z_e_mean: N x N x d_lat  (we use mean of edge latent for decoder input)
% E:        N x N  original edge features (carried through for attention)

% Summarise edge latent to scalar per (i,j): mean over latent dim
if ndims(z_e_mean) == 3
    E_lat_scalar = mean(z_e_mean, 3);  % N x N
else
    E_lat_scalar = z_e_mean;
end
% Combine original correlation with latent edge signal
E_dec = E + 0.1 * E_lat_scalar;

H = z_n;
for bi = 1:n_dec
    H = transConvBlock(p.dec(bi), H, E_dec, nh);
end

% 2-layer MLP to project back to d_emb for GRN inference
Z_prime = leakyrelu(H * p.dec_W1 + p.dec_b1, 0.01);
Z_prime = Z_prime * p.dec_W2 + p.dec_b2;
E_out = E_dec;
end

% -----------------------------------------------------------------------
function Y = transConvBlock(bp, X, E, nh)
% Transformer Convolution block with edge-augmented multi-head attention
% X:  N x d    node features
% E:  N x N    scalar edge weights (Pearson correlation)
% Y:  N x d    updated node features
%
% Attention:  a_ij = softmax( (q_i'(k_j + e_ij*We_k)) / sqrt(d_h) )
% Aggregate:  h_i  = sum_j( a_ij * (v_j + e_ij*We_v) )

[~, d] = size(X);
d_h    = d / nh;

Q = X * bp.Wq;    % N x d
K = X * bp.Wk;
V = X * bp.Wv;

% Compute edge contribution to keys: r_i = q_i^T We_k (scalar per node)
r = Q * bp.We_k';                        % N x 1 (r_i = <q_i, We_k>)

% Attention scores (N x N): S_ij = (q_i^T k_j + e_ij * r_i) / sqrt(d_h)
S = (Q * K' + E .* r) / sqrt(d_h);      % N x N
A = rowsoftmax(S);                       % row-wise (N x N)

% Value aggregation with edge contributions
E_agg = sum(A .* E, 2);                  % N x 1: sum_j(a_ij * e_ij)
H = A * V + E_agg .* bp.We_v;           % N x d

% Output projection + residual + layer norm 1
H  = H * bp.Wo;
X1 = layernorm(X + H, bp.ln1g, bp.ln1b);

% Feed-forward network
ff = leakyrelu(X1 * bp.ffW1 + bp.ffb1, 0.01);
ff = ff * bp.ffW2 + bp.ffb2;

% Residual + layer norm 2
Y = layernorm(X1 + ff, bp.ln2g, bp.ln2b);
end

% -----------------------------------------------------------------------
function probs = grnInfer(Z_prime, E_out)
% Inner product decoder for edge probabilities
% Z_prime: N x d_emb  decoded node representations
% probs:   N x N      probability matrix (after sigmoid)
logits = Z_prime * Z_prime';            % N x N inner product
% Add pooled edge feature (mean of edge attention across dim)
logits = logits + mean(E_out, 3);       % broadcast if E_out is 3-D
probs  = sigmoid(logits);
end

% -----------------------------------------------------------------------
function A = rowsoftmax(S)
% Numerically stable row-wise softmax (works with dlarray)
Sm = S - max(S, [], 2);
E  = exp(Sm);
A  = E ./ sum(E, 2);
end

% -----------------------------------------------------------------------
function LN = layernorm(X, g, b)
% Layer normalization along feature dimension (last/column)
% X: N x d,  g,b: 1 x d
mu = mean(X,2);
sg = std(X,0,2);  sg(sg<1e-8) = 1;
LN = ((X - mu) ./ sg) .* g + b;
end


% ======================================================================
%  VAE UTILITIES
% ======================================================================

function z = reparameterise(mu, lv)
% Gaussian reparameterisation trick: z = mu + exp(lv/2) * eps
eps = randn(size(mu), 'like', extractdata(mu));
z   = mu + exp(lv ./ 2) .* dlarray(eps);
end

% ======================================================================
%  LOSS
% ======================================================================

function [loss_total, loss_bce, loss_kl] = computeLoss(probs, y_true, mu_n, lv_n, beta)
% BCE reconstruction + normalised KL divergence
% probs, y_true: N x N   (predicted probs and binary labels)
% mu_n, lv_n:   N x d_lat (node latent parameters)
% beta: KL weight (1 / N_sub recommended)

EPS = 1e-7;
p  = max(min(probs, 1-EPS), EPS);     % clamp for stability
y  = dlarray(y_true);

% Dynamic negative sampling: match number of negatives to positives
n_pos = nnz(y_true > 0.5);
n_neg = nnz(y_true < 0.5);
if n_pos > 0 && n_neg > 0
    % Balanced binary cross-entropy via masked means
    pos_mask = logical(y_true > 0.5);
    neg_mask = ~pos_mask;
    bce_pos = -mean(log(p(pos_mask)));
    bce_neg = -mean(log(1 - p(neg_mask)));
    loss_bce = 0.5 * (bce_pos + bce_neg);
else
    loss_bce = -mean(y .* log(p) + (1-y) .* log(1-p), 'all');
end

% KL divergence (per-node, normalised)
loss_kl = -0.5 * mean( 1 + lv_n - mu_n.^2 - exp(lv_n), 'all');

loss_total = loss_bce + beta * loss_kl;
end


% ======================================================================
%  TRAINING
% ======================================================================

function params = trainModel(params, node_feats, sg_edges, sg_labels, ...
        N_sub, d_emb, d_lat, nh, n_enc, n_dec, ...
        n_epochs, batch_sz, lr, l1w, verbose)

n_sg   = numel(node_feats);
beta   = 1 / N_sub;
avg_g  = [];  avg_sq = [];   % Adam moment estimates
iter   = 1;

if verbose
    fprintf('[GRNFormer] Training: %d epochs, batch=%d, lr=%.4f\n', n_epochs, batch_sz, lr);
end

for ep = 1:n_epochs
    perm     = randperm(n_sg);
    ep_loss  = 0;
    n_batches = 0;

    for b_start = 1:batch_sz:n_sg
        b_idx = perm(b_start : min(b_start+batch_sz-1, n_sg));

        batch_nf  = node_feats(b_idx);
        batch_eg  = sg_edges(b_idx);
        batch_lb  = sg_labels(b_idx);

        % Compute loss + gradients via automatic differentiation
        [loss_val, grads] = dlfeval(@batchLossFn, params, ...
            batch_nf, batch_eg, batch_lb, d_emb, d_lat, nh, n_enc, n_dec, beta, l1w);

        [params, avg_g, avg_sq] = adamupdate(params, grads, avg_g, avg_sq, iter, lr);
        iter = iter + 1;

        ep_loss   = ep_loss + extractdata(loss_val);
        n_batches = n_batches + 1;
    end

    if verbose && (mod(ep,10)==0 || ep==1)
        fprintf('  epoch %3d/%d  loss=%.5f\n', ep, n_epochs, ep_loss/max(n_batches,1));
    end
end
end

% -----------------------------------------------------------------------
function [loss, grads] = batchLossFn(params, batch_nf, batch_eg, batch_lb, ...
                                      d_emb, ~, nh, n_enc, n_dec, beta, l1w)
% Called inside dlfeval -- computes batch loss and gradients
total_loss = dlarray(0);

for k = 1:numel(batch_nf)
    nf  = dlarray(batch_nf{k});        % N x d_in
    E   = dlarray(batch_eg{k});        % N x N
    lbl = batch_lb{k};                 % N x N  (double, constant)

    % --- Gene-Transcoder ---
    Z_emb = geneTranscoderDl(params.gt, nf, d_emb, nh);

    % --- GraViTAE encoder ---
    [mu_n, lv_n, mu_e, ~] = gravitaeEncode(params, Z_emb, E, n_enc, nh);

    % --- Reparameterise ---
    z_n = reparameterise(mu_n, lv_n);

    % --- GraViTAE decoder ---
    [Z_prime, E_out] = gravitaeDecode(params, z_n, mu_e, E, n_dec, nh);

    % --- GRN inference ---
    probs = grnInfer(Z_prime, E_out);

    % --- Loss ---
    [lk, ~, ~] = computeLoss(probs, lbl, mu_n, lv_n, beta);
    total_loss  = total_loss + lk;
end

total_loss = total_loss / numel(batch_nf);

% L1 regularisation on all weight matrices
l1_reg = dlarray(0);
fnames = fieldnames(params.gt);
for f = 1:numel(fnames)
    w = params.gt.(fnames{f});
    l1_reg = l1_reg + sum(abs(w),'all');
end
for bi = 1:numel(params.enc)
    fn = fieldnames(params.enc(bi));
    for f = 1:numel(fn)
        w = params.enc(bi).(fn{f});
        l1_reg = l1_reg + sum(abs(w),'all');
    end
end

total_loss = total_loss + l1w * l1_reg;

grads = dlgradient(total_loss, params);
loss  = total_loss;
end

% -----------------------------------------------------------------------
function Z_emb = geneTranscoderDl(p, nf, d_emb, nh)
% dlarray-compatible Gene-Transcoder (used inside dlfeval)
[N, T_cell] = size(nf);

% Pointwise Conv1D: project each cell's scalar → d_emb
nf_flat = reshape(nf, N*T_cell, 1);
tok_flat = nf_flat * p.W_conv + p.b_conv;        % (N*T) x d_emb
tokens   = reshape(tok_flat, N, T_cell, d_emb);  % N x T x d_emb

% Self-attention over cells per gene (all genes processed jointly)
% Flatten N and T: treat as (N*T)-length sequence split into N blocks of T
d_h = d_emb / nh;
embeds = dlarray(zeros(N, d_emb));

% Process gene by gene (inner loop over N=100 is acceptable)
for g = 1:N
    tok_g = squeeze(tokens(g,:,:));              % T_cell x d_emb
    Q = tok_g * p.Wq;
    K = tok_g * p.Wk;
    V = tok_g * p.Wv;
    A = rowsoftmax((Q * K') / sqrt(d_h));
    H = A * V * p.Wo;
    H1 = layernorm(tok_g + H, p.ln1g, p.ln1b);
    ff = leakyrelu(H1 * p.ffW1 + p.ffb1, 0.01);
    ff = ff * p.ffW2 + p.ffb2;
    out_g = layernorm(H1 + ff, p.ln2g, p.ln2b);
    embeds(g,:) = mean(out_g, 1);
end
Z_emb = embeds;
end


% ======================================================================
%  INFERENCE (no reparameterisation noise, returns double probs)
% ======================================================================

function probs = forwardInfer(params, nf, E, d_emb, ~, nh, n_enc, n_dec)
% Run full forward pass without gradient tracking (extractdata at end)
nf_dl = dlarray(nf);
E_dl  = dlarray(E);

Z_emb = geneTranscoder(params.gt, nf_dl, d_emb, nh);

[mu_n, ~, mu_e, ~] = gravitaeEncode(params, Z_emb, E_dl, n_enc, nh);

% Use mean (no noise) during inference
[Z_prime, E_out] = gravitaeDecode(params, mu_n, mu_e, E_dl, n_dec, nh);

prob_dl = grnInfer(Z_prime, E_out);
probs   = extractdata(prob_dl);
end

function Q = mpsqmi(mps, d)
% MPSQMI  Pairwise mutual information from an MPS Born machine's marginals.
%
%   Q = qtm.mpsqmi(mps, d)
%
%   Returns the N x N pairwise mutual information (bits) of the classical
%   distribution P(s) = |psi(s)|^2 / Z encoded by the trained MPS: it forms
%   the one- and two-site marginals of P and takes their Shannon MI -- the
%   same quantity qtm.BinPairMI estimates from data, here read from the model
%   instead.
%
%   NOTE: despite the name, this is CLASSICAL mutual information, not von
%   Neumann quantum mutual information. The Born state |psi> = sum_s sqrt(P(s))|s>
%   is a coherent superposition, not the classical mixture sum_s P(s)|s><s|;
%   its reduced density matrices carry off-diagonal coherences, so the true
%   QMI (S(rho_i) + S(rho_j) - S(rho_ij)) differs from the Shannon MI of the
%   marginals -- e.g. for perfectly correlated genes the Bell state gives 2
%   bits of QMI vs 1 bit of classical MI. Summing P over the other sites here
%   discards those coherences by construction. Computing true QMI requires
%   building rho_i, rho_ij from the sqrt(P) amplitudes and diagonalizing;
%   see +qtm/tngrn_quantum_notes.md and the RDM-route TODO below.
%
%   Panel-scale exact readout: contracts the full d^N joint. For larger panels
%   replace the joint contraction with a canonical-form reduced-density-matrix
%   route (bring the orthogonality centre between i and j for adjacent sites;
%   contract the intervening transfer matrices for distant ones).
%
%   Interventional note: if the MPS was clamped at a site (virtual knockout via
%   ten/i_mps_clamp), the amplitudes are unnormalised P(s, g=code); the
%   P = P./sum(P) step below renormalises to the conditional P(s | g=code), so
%   Q is then the post-knockout MI with no extra bookkeeping.
%
%   See also qtm.BinPairMI, qtm.MI_construction, ten.tngrn.

arguments
    mps cell
    d (1,1) double = 2
end

N    = numel(mps);
grid = i_allconfigs(N, d);          % [d^N x N], codes 1..d
psi  = i_amp(mps, grid);            % amplitudes
P    = psi .^ 2;
P    = P ./ sum(P);                 % normalise (also renormalises a clamped MPS)

ent = @(pv) -sum(pv(pv > 0) .* log2(pv(pv > 0)));
Q   = zeros(N);
for i = 1:N
    Pi = accumarray(grid(:, i), P, [d 1]);
    for j = i + 1:N
        Pj  = accumarray(grid(:, j), P, [d 1]);
        Pij = accumarray(grid(:, [i j]), P, [d d]);
        Q(i, j) = ent(Pi) + ent(Pj) - ent(Pij(:));
        Q(j, i) = Q(i, j);
    end
end
end % mpsqmi


% ------------------------------------------------------------------------
function grid = i_allconfigs(N, d)
grid = zeros(d ^ N, N);
for k = 1:N
    grid(:, k) = repmat(repelem((1:d).', d ^ (N - k)), d ^ (k - 1), 1);
end
end

function psi = i_amp(mps, cfg)
% Batched amplitude of configs [B x N] (codes 1..d) by left-to-right contract.
B = size(cfg, 1);
L = ones(B, 1);
for k = 1:numel(mps)
    Ak = mps{k}; [Dl, ~, Dr] = size(Ak);
    Ln = zeros(B, Dr);
    for s = 1:size(Ak, 2)
        m = cfg(:, k) == s; if ~any(m), continue; end
        Ln(m, :) = L(m, :) * reshape(Ak(:, s, :), [Dl, Dr]);
    end
    L = Ln;
end
psi = L(:, 1);
end

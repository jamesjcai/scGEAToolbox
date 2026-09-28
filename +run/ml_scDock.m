function ok = ml_scDock()
%ML_SCDOCK  Put the bundled scDock functions on the MATLAB path.
%   OK = RUN.ML_SCDOCK() adds external/ml_scDock to the path so that
%   SC_DOCK_CCC, SC_DOCK_VINA, SC_DOCK_GENE2PDB and SC_DOCK_RNA can be called.
%   Returns true when the folder was found.
%
%   GLY.WEIGHT and GLY.SHIELD are documented against SC_DOCK_CCC and
%   SC_DOCK_VINA, but those live in external/ rather than in the toolbox proper
%   and are not on the path by default. Call this first:
%
%     run.ml_scDock();
%     ccc = sc_dock_ccc(X, g, ctype);
%     T   = gly.weight(ccc.T_interactions, X, g, ctype);
%
%   Follows the same convention as RUN.ML_COMBAT: the path is left in place for
%   the rest of the session, and nothing is added in deployed builds.
%
% see also: SC_DOCK_CCC, GLY.WEIGHT, GLY.SHIELD, RUN.ML_COMBAT

pw1 = fileparts(mfilename('fullpath'));
target = fullfile(pw1, '..', 'external', 'ml_scDock');

ok = isfolder(target);
if ~ok
    warning("RUN:ML_SCDOCK:NotFound", ...
        "scDock folder not found at %s. The sc_dock_* functions are unavailable.", ...
        target);
    return;
end

if ~(ismcc || isdeployed)
    addpath(target);
end

end

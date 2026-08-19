function cfg = i_xctcfg(res, tag)
%I_XCTCFG  Assemble the TEN.I_XCTCORE configuration from parsed arguments.
%   CFG = TEN.I_XCTCFG(RES, TAG) packs the shared name-value results RES, as
%   returned by the inputParser of TEN.SCTENIFOLDXCT or
%   TEN.SCTENIFOLDXCT_GLYCO, into the struct TEN.I_XCTCORE expects. TAG is the
%   string used to prefix progress messages.
%
%   The correspondence hook and provenance flag default to the published
%   behaviour: weight-1 cognate pairs and no extra output columns. Callers that
%   want otherwise overwrite cfg.corrfcn and cfg.provenance after this returns.
%
%   Kept separate so the two entry points cannot disagree about defaults.
%
% see also: TEN.I_XCTCORE, TEN.SCTENIFOLDXCT, TEN.SCTENIFOLDXCT_GLYCO

arguments
    res (1, 1) struct
    tag (1, 1) string
end

cfg = struct();
cfg.n_dim = round(res.n_dim);
cfg.mu = res.mu;
cfg.pval = res.pval;
cfg.verbose = res.verbose;
cfg.w12mode = string(res.w12mode);
cfg.alpha = res.alpha;
cfg.grnoffset = res.grnoffset;
cfg.useparallel = logical(res.useparallel);

cfg.slv = struct( ...
    'name',    string(res.solver), ...
    'n_steps', res.n_steps, ...
    'lr',      res.lr, ...
    'seed',    res.seed);

cfg.tag = tag;
cfg.corrfcn = [];
cfg.provenance = false;

end % i_xctcfg

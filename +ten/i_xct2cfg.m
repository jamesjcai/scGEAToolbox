function cfg = i_xct2cfg(res, tag)
%I_XCT2CFG  Assemble the TEN.I_XCT2CORE configuration from parsed arguments.
%   CFG = TEN.I_XCT2CFG(RES, TAG) packs the shared name-value results RES, as
%   returned by the inputParser of TEN.SCTENIFOLDXCT2 or
%   TEN.SCTENIFOLDXCT2_GLYCO, into the struct TEN.I_XCT2CORE expects. TAG is the
%   string used to prefix progress messages.
%
%   The correspondence hook, provenance flag, precomputed-GRN slots and
%   candidate scope default to the published behaviour: weight-1 cognate
%   pairs in both samples, no extra output columns, a freshly-built network
%   per cell type per sample, and an output restricted to the L-R database
%   (candidates="database"). Callers that want otherwise overwrite
%   cfg.corrfcn, cfg.provenance, cfg.grn_s1, cfg.grn_t1, cfg.grn_s2, cfg.grn_t2,
%   cfg.candidates and cfg.topN after this returns.
%
%   Kept separate so the two entry points cannot disagree about defaults. The
%   differential path has no solver choice and no p-value cutoff, so this is not
%   TEN.I_XCTCFG with a different tag.
%
% see also: TEN.I_XCT2CORE, TEN.SCTENIFOLDXCT2, TEN.SCTENIFOLDXCT2_GLYCO,
%           TEN.I_XCTCFG

arguments
    res (1, 1) struct
    tag (1, 1) string
end

cfg = struct();
cfg.n_dim = round(res.n_dim);
cfg.mu = res.mu;
cfg.verbose = res.verbose;
cfg.w12mode = string(res.w12mode);
cfg.alpha = res.alpha;
cfg.grnoffset = res.grnoffset;
cfg.useparallel = logical(res.useparallel);

cfg.tag = tag;
cfg.corrfcn = [];
cfg.provenance = false;
cfg.grn_s1 = [];
cfg.grn_t1 = [];
cfg.grn_s2 = [];
cfg.grn_t2 = [];
cfg.grn_s1_processed = false;
cfg.grn_t1_processed = false;
cfg.grn_s2_processed = false;
cfg.grn_t2_processed = false;
cfg.candidates = "database";
cfg.topN = 200;

end % i_xct2cfg

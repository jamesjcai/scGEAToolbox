function cfg = i_xct2cfg(res, tag)
%I_XCT2CFG  Assemble the TEN.I_XCT2CORE configuration from parsed arguments.
%   CFG = TEN.I_XCT2CFG(RES, TAG) packs the shared name-value results RES, as
%   returned by the inputParser of TEN.SCTENIFOLDXCT2 or
%   TEN.SCTENIFOLDXCT2_GLYCO, into the struct TEN.I_XCT2CORE expects. TAG is the
%   string used to prefix progress messages.
%
%   The correspondence hook and provenance flag default to the published
%   behaviour: weight-1 cognate pairs in both samples and no extra output
%   columns. Callers that want otherwise overwrite cfg.corrfcn and
%   cfg.provenance after this returns.
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

end % i_xct2cfg

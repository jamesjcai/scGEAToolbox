function T = lectinmap(mapType)
%LECTINMAP  Glyco-module to lectin / LR-pair cognate maps.
%   T = GLY.LECTINMAP() returns the Channel-A map linking each glycan
%   biosynthesis module (the "ligand" side, whose per-cell activity is scored
%   by GLY.STATE) to the cognate lectin module (the "receptor" side, whose
%   member genes are read from GLY.GENESETS). It drives the glyco-lectin
%   communication channel in GLY.CCC.
%
%   T = GLY.LECTINMAP("lr") returns the Channel-B map linking canonical
%   ligand-receptor families (matched by symbol pattern) to the glyco module
%   that modulates their signalling and the compartment (sender/receiver) the
%   modulation acts on. It drives the glyco-weighting layer in GLY.WEIGHT.
%
%   All module names refer to sets defined in GLY.GENESETS.
%
%   OUTPUTS (mapType = "lectin", default):
%     T - table with columns:
%         Epitope      - short name of the glycan determinant
%         GlycoModules - comma-separated biosynthesis module names (sender);
%                        multiple modules are combined by their mean per-cell
%                        score to form the ligand signal
%         LectinModule - cognate lectin/reader module name (receiver)
%         Description  - human-readable description
%
%   OUTPUTS (mapType = "lr"):
%     T - table with columns:
%         LigandPattern   - case-insensitive regexp matched against the ligand
%         ReceptorPattern - case-insensitive regexp matched against the receptor
%         GlycoModule     - modulating glyco module name
%         Compartment     - "sender" or "receiver": which cell's glyco score
%                           weights the edge
%         Exponent        - +1 where the glycan ENABLES the interaction, -1
%                           where it MASKS it
%         Description     - human-readable description
%
%     Several rows may match the same ligand-receptor pair, and all of them
%     apply: GLY.LRWEIGHT multiplies their factors. This replaced a
%     first-matching-rule-wins scheme under which a pair could carry only one
%     module, which made a per-pair dependency fingerprint inexpressible.
%
% see also: GLY.CCC, GLY.WEIGHT, GLY.GENESETS, GLY.STATE

arguments
    mapType (1, 1) string = "lectin"
end

switch lower(mapType)
    case {"lectin", "a", "channela"}
        T = i_channelA();
    case {"lr", "b", "channelb"}
        T = i_channelB();
    otherwise
        error("LECTINMAP:BadType", ...
            "mapType must be ""lectin"" (Channel A) or ""lr"" (Channel B).");
end

end


%% ---- Channel A: glyco biosynthesis module -> cognate lectin module ----
function T = i_channelA()
% Each row: {Epitope, GlycoModules(sender), LectinModule(receiver), Description}.
spec = {

"sialic_acid", ...
"Glyco_sialylation", ...
"Glyco_siglecs", ...
"Sialylated glycans produced by sialyltransferases are read by Siglecs";

"ganglioside_sialic_acid", ...
"Glyco_glycosphingolipid_biosynthesis,Glyco_sialylation", ...
"Glyco_siglecs", ...
"Sialylated glycosphingolipids (gangliosides) engage Siglecs";

"sialyl_LewisX", ...
"Glyco_sialylation,Glyco_fucosylation", ...
"Glyco_selectins_and_ligands", ...
"Sialyl-Lewis-x (sialylation + fucosylation) is the selectin ligand";

"galactoside_polyLacNAc", ...
"Glyco_galactosylation,Glyco_N_glycan_processing_Golgi", ...
"Glyco_galectins", ...
"Beta-galactosides / poly-LacNAc on branched N-glycans are galectin ligands";

"fucose_mannose", ...
"Glyco_fucosylation", ...
"Glyco_C_type_lectins", ...
"Fucosylated / high-mannose glycans are recognised by C-type lectins";

};

T = cell2table(spec, VariableNames=["Epitope", "GlycoModules", ...
    "LectinModule", "Description"]);
T.Epitope      = string(T.Epitope);
T.GlycoModules = string(T.GlycoModules);
T.LectinModule = string(T.LectinModule);
T.Description  = string(T.Description);
end


%% ---- Channel B: LR family -> modulating glyco module(s) & compartment ----
function T = i_channelB()
% Each row: {LigandPattern, ReceptorPattern, GlycoModule, Compartment,
%            Exponent, Description}.
%
% EVERY matching row applies, multiplicatively - see GLY.LRWEIGHT. A pair
% may therefore carry several modules, which is the point: sialyl-Lewis-x needs
% sialylation AND fucosylation, and EGFR surface residency needs branching AND
% the galectin lattice that binds the branches. Before this the first matching
% rule won and every later one was silently discarded.
%
% EXPONENT is the sign of the dependence, following the a_ijk of the
% perturbation formulation: +1 where the glycan ENABLES the interaction, -1
% where it MASKS it. Polysialylation of NCAM is the worked example of the
% negative case - it attenuates NCAM-FGFR signalling rather than supporting it.
% Do not add a -1 row without a mechanism; a sign error inverts the prediction
% rather than weakening it.
spec = {

"^FGF\d", "^FGFR", ...
"Glyco_heparan_sulfate_biosynthesis", "receiver", 1, ...
"Heparan sulfate is an obligate co-receptor for FGF-FGFR signalling";

"^WNT", "^FZD", ...
"Glyco_heparan_sulfate_biosynthesis", "receiver", 1, ...
"Heparan sulfate proteoglycans shape Wnt-Frizzled gradients and binding";

"^(CXCL|CCL)", "^(CXCR|CCR|ACKR)", ...
"Glyco_heparan_sulfate_biosynthesis", "receiver", 1, ...
"Chemokines are presented and stabilised by cell-surface heparan sulfate";

"^(CXCL|CCL)", "^(CXCR|CCR|ACKR)", ...
"Glyco_carbohydrate_sulfotransferases", "receiver", 1, ...
"Chemokine binding to heparan sulfate is sulfation-pattern dependent";

"^(PTN|MDK)", "^SDC", ...
"Glyco_heparan_sulfate_biosynthesis", "receiver", 1, ...
"Pleiotrophin and midkine bind the heparan sulfate chains of syndecans";

"^(PTN|MDK)", "^PTPRZ1", ...
"Glyco_chondroitin_dermatan_sulfate", "receiver", 1, ...
"PTPRZ1 carries chondroitin sulfate chains that PTN and MDK bind";

"^(DLL|JAG)", "^NOTCH", ...
"Glyco_Notch_EGF_O_glycosylation", "receiver", 1, ...
"O-fucose/Fringe glycosylation of Notch EGF repeats tunes ligand selectivity";

"^(EGF|TGFA|AREG|EREG|HBEGF)", "^(EGFR|ERBB)", ...
"Glyco_N_glycan_processing_Golgi", "receiver", 1, ...
"N-glycan branching (MGAT5) controls EGFR-family surface residency";

"^(EGF|TGFA|AREG|EREG|HBEGF)", "^(EGFR|ERBB)", ...
"Glyco_galectins", "receiver", 1, ...
"Galectins cross-link the branched N-glycans into the lattice that retains EGFR";

"^TGFB", "^TGFBR", ...
"Glyco_N_glycan_processing_Golgi", "receiver", 1, ...
"N-glycan branching modulates TGF-beta receptor surface retention";

"^(SELPLG|CD34|PODXL|PODXL2|GLG1)", "^SEL[EPL]", ...
"Glyco_sialylation", "sender", 1, ...
"Selectin ligands require the sialic acid of sialyl-Lewis-x on the presenting cell";

"^(SELPLG|CD34|PODXL|PODXL2|GLG1)", "^SEL[EPL]", ...
"Glyco_fucosylation", "sender", 1, ...
"Selectin ligands equally require the fucose of sialyl-Lewis-x (FUT7)";

"^NCAM", "^FGFR", ...
"Glyco_sialylation", "sender", -1, ...
"Polysialylation of NCAM (ST8SIA2/4) MASKS NCAM-FGFR signalling";

};

T = cell2table(spec, VariableNames=["LigandPattern", "ReceptorPattern", ...
    "GlycoModule", "Compartment", "Exponent", "Description"]);
T.LigandPattern   = string(T.LigandPattern);
T.ReceptorPattern = string(T.ReceptorPattern);
T.GlycoModule     = string(T.GlycoModule);
T.Compartment     = string(T.Compartment);
T.Exponent        = double(cell2mat(num2cell(T.Exponent)));
T.Description     = string(T.Description);
end

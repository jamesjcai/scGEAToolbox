function [paragraph, refkeys] = i_damethodstext(T, options)
%I_DAMETHODSTEXT Methods text for a differential cell-type abundance run.
%
%   [paragraph, refkeys] = pkg.i_damethodstext(T, Name=Value)
%
%   Describes what SC_DIFFABUNDANCE did, reading the test it ran from
%   T.Properties.Description, so the text cannot name a different test
%   from the one in the table.
%
%   NAME-VALUE ARGUMENTS
%     CellTypeVar   name of the cell-type variable
%     GroupVar      name of the condition variable
%     SampleVar     name of the replicate variable, "" for none
%     Groups        the conditions compared, in the order tested
%     NumSamples    number of distinct replicate samples
%
%   See also SC_DIFFABUNDANCE, PKG.I_METHODSREPORT.

arguments
    T table
    options.CellTypeVar (1,1) string = "cell type"
    options.GroupVar (1,1) string = "condition"
    options.SampleVar (1,1) string = ""
    options.Groups (1,:) string = strings(1, 0)
    options.NumSamples (1,1) double = NaN
end

desc = string(T.Properties.Description);
grpText = "conditions";
if ~isempty(options.Groups)
    grpText = strjoin("""" + options.Groups + """", ", ");
end
lead = sprintf("Differences in the proportions of %d %s groups between " + ...
    "%s (%s) were tested with sc_diffabundance in scGEAToolbox.", ...
    height(T), lower(options.CellTypeVar), grpText, options.GroupVar);

refkeys = ["benjamini1995", "cai2020"];
if contains(desc, "propeller")
    if contains(desc, "asin")
        tr = "an arcsine square-root transformation";
    else
        tr = "a logit transformation (0.5 added to every count)";
    end
    if contains(desc, "moderated F")
        test = "an empirical Bayes moderated F-test";
    else
        test = "an empirical Bayes moderated t-test";
    end
    body = sprintf("Following propeller, cell-type proportions were " + ...
        "computed per replicate sample (%s) and condition, variance-" + ...
        "stabilised with %s, and fitted per cell type with a linear " + ...
        "model on condition", options.SampleVar, tr);
    if contains(desc, "blocked by sample")
        body = body + " with sample as a blocking factor";
    end
    body = body + "; the condition effect was tested with " + test + ...
        " that shares variance information across cell types.";
    if ~isnan(options.NumSamples)
        body = body + sprintf(" %d samples were analysed.", options.NumSamples);
    end
    refkeys = ["phipson2022", "smyth2004", refkeys];
elseif contains(desc, "permutation")
    body = "Each cell was treated as an observation and the difference " + ...
        "in proportion was compared with that under randomly permuted " + ...
        "condition labels; this does not account for sample-to-sample " + ...
        "variation.";
elseif contains(desc, "Fisher")
    body = "Each cell type was tested with Fisher's exact test on its cell " + ...
        "counts against all other cells; this does not account for " + ...
        "sample-to-sample variation.";
else
    body = sprintf("Per-sample proportions (%s) were compared with a " + ...
        "%s.", options.SampleVar, extractAfter(desc, ": "));
end
paragraph = lead + " " + body + " P-values were adjusted across cell " + ...
    "types by the Benjamini-Hochberg procedure, and changes with FDR < " + ...
    "0.05 were reported.";
end

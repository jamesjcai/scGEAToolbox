function [S] = e_doubletcallsummary(isDoublet, info)
%E_DOUBLETCALLSUMMARY What a Scrublet run can honestly be reported as.
%   S = pkg.e_doubletcallsummary(isDoublet, info) turns SC_SCRUBLET's
%   outputs into the message a user should see, so that the three states
%   it can end in are told apart instead of collapsed.
%
%   S.Status        "noThreshold" | "none" | "weak" | "called"
%   S.Title         dialog title
%   S.Message       what to say
%   S.IsWarning     true when the run must not be read as a clean result
%   S.OfferRemoval  true when offering to delete the called cells makes
%                   sense
%
%   INFO may be empty, as it is for the Python backend, in which case
%   only the calls themselves are available.
%
%   Three things this exists to stop GUI.CALLBACK_DOUBLETDETECTION doing.
%
%   SC_SCRUBLET returns an all-false ISDOUBLET for two unrelated
%   reasons: a threshold was found and no cell cleared it, which is a
%   negative result; or the simulated scores were not bimodal, no
%   threshold could be found and nothing was classified at all, which it
%   signals with INFO.THRESHOLD = NaN and a console warning. Both were
%   reported as "No doublet detected.", so a failed screen read as a
%   clean QC pass.
%
%   SC_SCRUBLET also warns, again only to the console, when fewer than
%   30% of the simulated doublets clear the threshold -- most of them
%   are then indistinguishable from single cells and the calls should
%   not be trusted. An App Designer user never sees the console.
%
%   And the reported "true rate" is detectedRate/detectableFraction,
%   which is only a rate while the denominator is large enough. On pure
%   Poisson noise, measured: 122 of 150 cells called (81.3%), 23% of
%   simulated doublets detectable, and the dialog announcing "The true
%   rate is therefore nearer 348.6%" before offering to delete those 122
%   cells.
%
%   See also SC_SCRUBLET, GUI.CALLBACK_DOUBLETDETECTION.

arguments
    isDoublet (:, 1) logical
    info = []
end

% SC_SCRUBLET's own cut-off for "the threshold rests on a weak split".
weakBelow = 0.3;

nCalled = sum(isDoublet);
S = struct("Status", "called", "Title", "Doublet Detection", ...
    "Message", "", "IsWarning", false, "OfferRemoval", nCalled > 0);

threshold = i_field(info, "threshold", NaN);
detectable = i_field(info, "detectableFraction", NaN);
detectedRate = i_field(info, "detectedDoubletRate", mean(isDoublet));

if ~isempty(info) && isnan(threshold)
    S.Status = "noThreshold";
    S.IsWarning = true;
    S.OfferRemoval = false;
    S.Message = ['No cell was tested. The simulated doublet scores ', ...
        'are not bimodal, so no threshold could be found and nothing ', ...
        'was classified -- this is not a negative result. It is what ', ...
        'a single homogeneous population looks like. Inspect the ', ...
        'simulated scores, or pass a threshold explicitly to ', ...
        'sc_scrublet.'];
    return
end

if nCalled == 0
    S.Status = "none";
    S.Message = 'No doublet detected.';
    return
end

pct = @(x) 100*x;
callLine = sprintf('%d cells called as doublets (%.1f%%).', ...
    nCalled, pct(detectedRate));

if isnan(detectable) || detectable <= 0
    % No INFO, or nothing to say about detectability.
    S.Message = callLine;
    return
end

detectLine = sprintf(['Only doublets of two different cell types are ', ...
    'detectable, and %.0f%% of the simulated ones were.'], pct(detectable));

overall = detectedRate/detectable;
if overall <= 1
    rateLine = sprintf('The true rate is therefore nearer %.1f%%.', ...
        pct(overall));
else
    % detectedRate/detectableFraction stops being a rate once it exceeds
    % 1. Printing it anyway produced "nearer 348.6%".
    rateLine = ['Correcting for that would put the rate above 100%, ', ...
        'which means the correction does not apply here: far more ', ...
        'cells were called than the method can actually resolve.'];
end

if detectable < weakBelow
    S.Status = "weak";
    S.IsWarning = true;
    S.Message = sprintf('%s\n\n%s\n\n%s\n\n%s', callLine, detectLine, ...
        rateLine, ['These calls should not be trusted. Most simulated ', ...
        'doublets are indistinguishable from single cells here, so ', ...
        'the threshold rests on a weak split rather than on a real ', ...
        'valley in the score histogram.']);
else
    S.Message = sprintf('%s\n\n%s %s', callLine, detectLine, rateLine);
end
end

function v = i_field(info, name, default)
if isstruct(info) && isscalar(info) && isfield(info, name) && ...
        isscalar(info.(name)) && isnumeric(info.(name))
    v = double(info.(name));
else
    v = default;
end
end

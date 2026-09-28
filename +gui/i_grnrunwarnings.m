function warnings = i_grnrunwarnings(method, x, numNetworks)
%I_GRNRUNWARNINGS The pre-run warnings that apply to one network run.
%   W = gui.i_grnrunwarnings(METHOD, X, NUMNETWORKS) returns a struct
%   array with fields Title, Message and Default, one element per warning
%   gui.i_confirmgrnrun should show, for METHOD (a row of
%   net.grnmethods()) on the transformed genes-by-cells matrix X:
%     - METHOD has a RawCountsNote and X still looks like raw counts
%       (whole-number values, i.e. "No" in the transform dialog);
%     - METHOD has a SecondsPerPair estimate and NUMNETWORKS networks
%       would take over a minute.
%   Kept apart from the dialogs so the rules can be tested without them.
%
%   See also gui.i_confirmgrnrun, net.grnmethods.

warnings = struct("Title", {}, "Message", {}, "Default", {});

if method.RawCountsNote ~= ""
    v = nonzeros(x);
    if all(v == round(v))
        warnings(end + 1) = struct("Title", method.Label + " input", ...
            "Message", "The expression matrix looks like raw counts. " + ...
            method.RawCountsNote, "Default", "Cancel");
    end
end

if ~isnan(method.SecondsPerPair)
    numGenes = size(x, 1);
    numPairs = numNetworks*numGenes*(numGenes - 1)/2;
    minutes = numPairs*method.SecondsPerPair/60;
    if minutes > 1
        warnings(end + 1) = struct("Title", upper(method.Key) + " run time", ...
            "Message", sprintf(['%s on %d genes scores %d gene pairs, ', ...
            'which takes about %.0f minutes (less with a parallel pool ', ...
            'open). Continue?'], upper(method.Key), numGenes, numPairs, minutes), ...
            "Default", "Continue");
    end
end
end

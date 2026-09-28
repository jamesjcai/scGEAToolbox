function s = i_enrichrstatus(errs)
%I_ENRICHRSTATUS Overall status of a batch of Enrichr calls.
%   s = i_enrichrstatus(errs) is 'completed' when no call failed, 'failed'
%   when every one did, and 'partial' otherwise, so one bad library call
%   does not hide the rest. ERRS holds each call's error message, '' for
%   success.
errs = string(errs);
if isempty(errs) || all(errs == "")
    s = 'completed';
elseif all(errs ~= "")
    s = 'failed';
else
    s = 'partial';
end
end

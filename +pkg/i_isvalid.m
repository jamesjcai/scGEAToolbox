function tf = i_isvalid(h)
%I_ISVALID Robust isvalid wrapper for compiled (standalone) MATLAB applications.
%   In deployed apps, isvalid can throw instead of returning false for
%   deleted or invalid handles. This wrapper catches that failure mode.
    if isempty(h)
        tf = false;
        return;
    end
    try
        tf = isvalid(h);
    catch
        tf = false;
    end
end

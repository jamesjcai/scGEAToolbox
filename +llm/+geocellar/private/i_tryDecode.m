function data = i_tryDecode(text)
data = {};
try
    parsed = jsondecode(text);
    if isstruct(parsed); parsed = num2cell(parsed); end
    if iscell(parsed) && ~isempty(parsed) && isstruct(parsed{1})
        data = parsed;
    end
catch
    % non-JSON input returns the default {} per this helper's contract
end
end

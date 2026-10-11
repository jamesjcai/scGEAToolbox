function s = i_thousands(x)
%I_THOUSANDS A count written with thousands separators, for prose.
%
%   s = pkg.i_thousands(12113) returns "12,113". A non-integer is returned
%   as STRING writes it. Used by the Methods-text builders, so every count
%   in a Methods section is written the same way.
%
%   See also PKG.I_METHODSREPORT.

if x ~= round(x)
    s = string(x);
    return
end
s = string(regexprep(sprintf('%d', round(x)), '(\d)(?=(\d{3})+$)', '$1,'));
end

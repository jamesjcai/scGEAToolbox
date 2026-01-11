function [szA, szB, szOut] = setupImplicitExpansion(szA, szB)
% Internal use only.

%   Copyright 2025 The MathWorks, Inc.

% Determine output size
ndOut = max(length(szA), length(szB));
szA(end+1:ndOut) = 1;
szB(end+1:ndOut) = 1;
szOut = zeros(1, ndOut);
for ii = 1:ndOut
    if szA(ii) == szB(ii)
        % Matched
        szOut(ii) = szA(ii);
    elseif szA(ii) == 1
        % Expand in B
        szOut(ii) = szB(ii);
    elseif szB(ii) == 1
        % Expand in A
        szOut(ii) = szA(ii);
    else
        error(message('MATLAB:sizeDimensionsMustMatch'))
    end
end
end
function code = buildGrayCode(numBits)
% Returns 2^numBits-by-numBits char array defining the Gray Code for
% the input number of bits. Two consecutive rows of code differ by 1 bit.

% Copyright 2023-2024 The MathWorks, Inc.
arguments
    numBits (1,1) {mustBeInteger, mustBePositive}
end

% Preallocate
code = repmat('0', 2^numBits, numBits);

code(:, 1) = repelem('01', 2^(numBits-1));
for ii=2:numBits
    block = repelem('0110', 2^(numBits-ii));
    code(:, ii) = repmat(block, 1, 2^(ii-2));
end
end
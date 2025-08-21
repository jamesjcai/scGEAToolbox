function bitStr = hex2bin(hexStr, numQubits)
% Internal use only. Convert hex string array to its binary representation
% with the specified number of qubits.

% This is only used for converting measurement samples returned from IBM.

% Copyright 2024 The MathWorks, Inc.

assert(isstring(hexStr) && isvector(hexStr))

hexStr = lower(hexStr);

assert(all(startsWith(hexStr, "0x")))

hexStr = extractAfter(hexStr, "0x");

maxNumHex = max(strlength(hexStr));
if maxNumHex < 13
    % Convert full hex to double-precision integer then to binary
    bitStr = string(dec2bin(hex2dec(hexStr), numQubits));
else
    % Convert each hex character to a length-4 binary string
    hexChar = char(pad(hexStr, maxNumHex, "left", "0"));
    [m,n] = size(hexChar);
    maxNumBits = 4*n;
    bitChar = reshape(dec2bin(hex2dec(hexChar(:)), 4), [m, n, 4]);
    bitChar = reshape(permute(bitChar, [1 3 2]), [m maxNumBits]);
    % Pad on the left or trim from the right so bitStr has numQubits
    if numQubits > maxNumBits
        bitChar = [repmat('0', [m numQubits-maxNumBits]) bitChar];
    elseif numQubits < maxNumBits
        bitChar = bitChar(:, end-numQubits+1:end);
    end
    bitStr = string(bitChar);
end
end




function gates = mcxgraycode(ctrls, trgt)
%

% Copyright 2022-2024 The MathWorks, Inc.

Nc = length(ctrls);
graycode = quantum.internal.gate.buildGrayCode(Nc);
gates = [];

% Construct a MCZ gate array by setting lambda so that R1^2^(Nc-1) = Z
lambda = pi/(2^(Nc-1));
last_code = graycode(1, :); 
for i = 2:length(graycode)
    code = graycode(i, :);

    % index of first set bit, i.e., left-most '1' in the code
    % Never empty, because only first row of graycode is all-zero
    ones_bits = strfind(code, '1');
    set_idx = ones_bits(1);

    % index of the changed bit
    % Always scalar, graycode is constructed so only one bit changes row-to-row.
    diff_idx = find(code ~= last_code);

    % flipped bit becomes the control
    if diff_idx ~= set_idx 
        gates = [gates; cxGate(ctrls(diff_idx), ctrls(set_idx))]; %#ok<*AGROW>
    % if equal, the next left-most '1' is used, if possible 
    elseif length(ones_bits)>=2
        next_idx = ones_bits(2);
        gates = [gates; cxGate(ctrls(next_idx), ctrls(set_idx))]; 
    end

    % check parity of the code 
    if (mod(count(code, '1'), 2) == 0)
        gates = [gates; cr1Gate(ctrls(set_idx), trgt, -lambda)]; 
    else
        gates = [gates; cr1Gate(ctrls(set_idx), trgt, lambda)]; 
    end
    last_code = code;
end

% Use the identity X = HZH to create the final MCX gate array
gates = [hGate(trgt); gates; hGate(trgt)];

end

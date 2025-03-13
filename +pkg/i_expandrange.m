function expandedGSMs = i_expandrange(gsmRange)
    % EXPANDGSMRANGE Expands a GSM range notation to individual GSMs
    %
    % Input:
    %   gsmRange - string in format 'GSMXXXXXXX-Y' where Y indicates 
    %              the last digit of the final GSM in the range
    %
    % Output:
    %   expandedGSMs - comma-separated string of individual GSMs
    
    % Parse the input
    parts = split(gsmRange, '-');
    
    if length(parts) ~= 2
        error('Input format should be GSMXXXXXXX-Y');
    end
    
    baseGSM = parts{1};
    endDigit = parts{2};
    
    % Extract the base GSM number
    if ~startsWith(baseGSM, 'GSM')
        error('The GSM should start with "GSM"');
    end
    
    basePrefix = 'GSM';
    baseNumberStr = baseGSM(4:end);
    baseNumber = str2double(baseNumberStr);
    
    % Get the final digit of the base GSM
    baseLastDigit = str2double(baseNumberStr(end));
    
    % Get the end digit
    endDigitValue = str2double(endDigit);
    
    if isnan(endDigitValue)
        error('The number after the hyphen should be a numeric digit');
    end
    
    % Calculate how many GSMs to generate
    % If the end digit is less than the base last digit, we assume it wraps around
    if endDigitValue >= baseLastDigit
        count = endDigitValue - baseLastDigit + 1;
    else
        count = endDigitValue + 10 - baseLastDigit + 1;
    end
    
    % Generate the sequence of GSMs
    gsmList = cell(count, 1);
    for i = 0:(count-1)
        gsmList{i+1} = sprintf('%s%d', basePrefix, baseNumber + i);
    end
    
    % Join the GSMs with commas
    expandedGSMs = strjoin(gsmList, ',');
end




%{
I need a MATLAB function to convert: 'GSM8122873-6' to 'GSM8122873,GSM8122874,GSM8122875,GSM8122876'
%}

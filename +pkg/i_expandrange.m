function expanded = i_expandrange(input)
    % EXPANDGSMSERIES Expands GSM series notation to comma-separated list
    % Examples:
    %   expandGSMSeries('GSM12323-6') returns 'GSM12323,GSM12324,GSM12325,GSM12326'
    %   expandGSMSeries('GSM12328-31') returns 'GSM12328,GSM12329,GSM12330,GSM12331'
    
    % Check if input contains a hyphen
    if ~contains(input, '-')
        expanded = input;
        return;
    end
    
    % Extract prefix and range
    [prefix, numStart, numEnd] = parseGSMRange(input);
    
    % Generate all numbers in the range
    numbers = numStart:numEnd;
    
    % Format each number with the prefix
    gsmList = cell(1, length(numbers));
    for i = 1:length(numbers)
        gsmList{i} = [prefix num2str(numbers(i))];
    end
    
    % Join with commas
    expanded = strjoin(gsmList, ',');
end

function [prefix, numStart, numEnd] = parseGSMRange(input)
    % Parse GSM range notation
    % Find position of hyphen
    hyphenPos = strfind(input, '-');
    
    % Extract parts before and after hyphen
    beforeHyphen = input(1:hyphenPos-1);
    afterHyphen = input(hyphenPos+1:end);
    
    % Extract numeric part from before hyphen
    prefix = '';
    numStartStr = '';
    for i = length(beforeHyphen):-1:1
        if isstrprop(beforeHyphen(i), 'digit')
            numStartStr = [beforeHyphen(i), numStartStr];
        else
            % First non-digit from right is the end of prefix
            prefix = beforeHyphen(1:i);
            break;
        end
    end
    
    % Convert starting number to integer
    numStart = str2double(numStartStr);
    
    % Determine end number based on how many digits are in the afterHyphen part
    numEndStr = afterHyphen;
    
    % If afterHyphen has fewer digits than numStartStr, use prefix digits
    if length(numEndStr) < length(numStartStr)
        prefixDigits = numStartStr(1:length(numStartStr)-length(numEndStr));
        numEndStr = [prefixDigits, numEndStr];
    end
    
    % Convert ending number to integer
    numEnd = str2double(numEndStr);
    
    % Sanity check - ensure numEnd is greater than numStart
    if numEnd < numStart
        error('Invalid range: end number is less than start number');
    end
end

%{
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





I need a MATLAB function to convert: 'GSM8122873-6' to 'GSM8122873,GSM8122874,GSM8122875,GSM8122876'

% In matlab, write a function to convert 'GSM12323-6' to 'GSM12323,GSM12324,GSM12325,GSM12326'. The function should also can convert 'GSM12328-31' to 'GSM12328,GSM12329,GSM12330,GSM12331'


%}

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


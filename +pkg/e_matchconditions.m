function mask = e_matchconditions(T, conditions, combine)
%E_MATCHCONDITIONS  Which rows of a table meet a set of conditions.
%
%   mask = pkg.e_matchconditions(T, conditions) returns a logical column,
%   true for the rows of table T that meet every condition. CONDITIONS is a
%   struct array with fields
%
%     Variable  name of a column of T
%     Operator  one of the operators below
%     Value     what the column is compared with: a number or text for a
%               numeric column (text is parsed with STR2DOUBLE), text for a
%               text column
%
%   mask = pkg.e_matchconditions(T, conditions, "any") is true for the rows
%   that meet at least one condition. "all", the default, means every one.
%   No conditions at all match every row.
%
%   Operators for numeric columns: "Equals", "Not equal", "Greater than",
%   "At least", "Less than", "At most". "Equals" allows for rounding, a
%   relative difference of 1e-9.
%   Operators for text columns: "Equals", "Not equal" (exact, case
%   sensitive), "Contains", "Starts with" (case insensitive), "Matches
%   regex" (REGEXP, case insensitive, anywhere in the text).
%   PKG.E_CONDITIONOPERATORS returns these lists.
%
%   A missing value (NaN, or missing text) meets no condition, not even
%   "Not equal".
%
%   See also PKG.E_CONDITIONOPERATORS, PKG.E_CELLCONDITIONTABLE,
%   GUI.CALLBACK_BRUSHCELLSBYCONDITIONS.

arguments
    T table
    conditions struct
    combine (1, 1) string {mustBeMember(combine, ["all", "any"])} = "all"
end

numericOperators = pkg.e_conditionoperators("numeric");
textOperators = pkg.e_conditionoperators("text");

numRows = height(T);
if isempty(conditions)
    mask = true(numRows, 1);
    return;
end

hits = false(numRows, numel(conditions));
for k = 1:numel(conditions)
    hits(:, k) = i_match(T, conditions(k), numericOperators, textOperators);
end
if combine == "all"
    mask = all(hits, 2);
else
    mask = any(hits, 2);
end
end

function hit = i_match(T, condition, numericOperators, textOperators)
name = string(condition.Variable);
if ~ismember(name, string(T.Properties.VariableNames))
    error("e_matchconditions:UnknownVariable", ...
        "The table has no column named ""%s"".", name);
end
x = T.(name);
op = string(condition.Operator);

if isnumeric(x) || islogical(x)
    if ~ismember(op, numericOperators)
        error("e_matchconditions:BadOperator", ...
            """%s"" cannot compare the numbers in ""%s"". Use one of: %s.", ...
            op, name, strjoin(numericOperators, ", "));
    end
    v = condition.Value;
    if ~isnumeric(v)
        v = str2double(string(v));
    end
    if ~isscalar(v) || isnan(v)
        error("e_matchconditions:BadValue", ...
            "Column ""%s"" holds numbers. Give a number to compare it with.", name);
    end
    x = double(x);
    switch op
        case "Equals"
            hit = abs(x - v) <= 1e-9*max(1, abs(v));
        case "Not equal"
            hit = abs(x - v) > 1e-9*max(1, abs(v));
        case "Greater than"
            hit = x > v;
        case "At least"
            hit = x >= v;
        case "Less than"
            hit = x < v;
        case "At most"
            hit = x <= v;
        otherwise
            % Every numeric operator is handled above.
            hit = false(size(x));
    end
    hit = hit & ~isnan(x);
else
    if ~ismember(op, textOperators)
        error("e_matchconditions:BadOperator", ...
            """%s"" cannot compare the text in ""%s"". Use one of: %s.", ...
            op, name, strjoin(textOperators, ", "));
    end
    x = string(x);
    v = string(condition.Value);
    switch op
        case "Equals"
            hit = x == v;
        case "Not equal"
            hit = x ~= v;
        case "Contains"
            hit = contains(x, v, IgnoreCase=true);
        case "Starts with"
            hit = startsWith(x, v, IgnoreCase=true);
        case "Matches regex"
            hit = ~cellfun(@isempty, regexpi(cellstr(x), char(v), "once"));
        otherwise
            % Every text operator is handled above.
            hit = false(size(x));
    end
    hit = hit & ~ismissing(x);
end
hit = hit(:);
end

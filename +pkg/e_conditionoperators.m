function ops = e_conditionoperators(kind)
%E_CONDITIONOPERATORS  The operators PKG.E_MATCHCONDITIONS accepts for a kind of column.
%
%   ops = pkg.e_conditionoperators("numeric") returns the operators for a
%   numeric column, and pkg.e_conditionoperators("text") those for a text
%   column, as string rows in the order a menu should list them.
%
%   See also PKG.E_MATCHCONDITIONS.

arguments
    kind (1, 1) string {mustBeMember(kind, ["numeric", "text"])}
end

if kind == "numeric"
    ops = ["Equals", "Not equal", "Greater than", "At least", "Less than", ...
        "At most"];
else
    ops = ["Equals", "Not equal", "Contains", "Starts with", "Matches regex"];
end
end

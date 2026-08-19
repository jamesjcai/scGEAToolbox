function p = i_normalizepath(p)
% I_NORMALIZEPATH Resolve "." and ".." in a path to a canonical absolute path
%
%   p = i_normalizepath(p) lexically normalizes P, collapsing "." and ".."
%   segments and returning an absolute path. It uses only core MATLAB (no
%   Java runtime), so it works on installations without a bundled JRE.
%
%   Example:
%   p = i_normalizepath("D:\GitHub\scGEAToolbox_dev\+run\..\external\script.py")
%   % -> "D:\GitHub\scGEAToolbox_dev\external\script.py"

if isstring(p) || ischar(p)
    p = char(p);
else
    error("Input must be a string or char array.");
end

if isempty(p)
    return;
end

% On Windows, treat both "/" and "\" as separators. On Unix, "\" is a valid
% filename character, so only "/" is a separator. Char literals are used
% deliberately: mixing a string scalar into strrep would promote P to a
% string and break the char indexing below.
if ispc
    p = strrep(p, '/', filesep);
end

% Make the path absolute (relative to the current folder) if needed.
if ~in_isabsolute(p)
    p = [pwd, filesep, p];
    if ispc
        p = strrep(p, '/', filesep);
    end
end

% Separate the root prefix so ".." never collapses past it:
%   Windows: "C:\", "C:", or UNC "\\server\share"
%   Unix:    leading "/"
if ispc
    prefix = regexp(p, '^(\\\\[^\\]+\\[^\\]+|[A-Za-z]:\\?)', "match", "once");
elseif startsWith(p, filesep)
    prefix = filesep;
else
    prefix = "";
end
rest = p(numel(prefix)+1:end);

% Collapse "." and ".." segments.
parts = strsplit(rest, filesep);
stack = cell(1, 0);
for k = 1:numel(parts)
    seg = parts{k};
    if isempty(seg) || strcmp(seg, ".")
        continue;
    elseif strcmp(seg, "..")
        if ~isempty(stack)
            stack(end) = [];
        end
    else
        stack{end+1} = seg; %#ok<AGROW>
    end
end

% Reassemble, guaranteeing exactly one separator after the root prefix.
prefix = char(prefix);
if ispc && ~isempty(prefix) && ~endsWith(prefix, filesep)
    prefix = [prefix, filesep]; % "C:" -> "C:\"
end
p = [prefix, strjoin(stack, filesep)];

% Drop a trailing separator except when the whole path is just the root.
if numel(p) > numel(prefix) && endsWith(p, filesep)
    p = p(1:end-1);
end
end

function tf = in_isabsolute(p)
if ispc
    tf = ~isempty(regexp(p, '^([A-Za-z]:\\|[A-Za-z]:/|\\\\)', "once"));
else
    tf = startsWith(p, "/");
end
end

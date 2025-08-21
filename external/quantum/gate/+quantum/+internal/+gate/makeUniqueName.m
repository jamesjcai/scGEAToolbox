function name = makeUniqueName(existing_names, varargin)
%

%   Copyright 2022 The MathWorks, Inc.

if isempty(varargin)
    basename = "cg";
else
    basename = varargin{1};
end

isUnique = 0;
id = 0;

while ~isUnique
   id = id+1;
   name = sprintf("%s%g", basename, id);
   isUnique = ~any(strcmp(existing_names, name));
end

end

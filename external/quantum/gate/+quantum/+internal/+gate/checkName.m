function checkName(name)
%

%   Copyright 2021-2022 The MathWorks, Inc.

if isequal(name, "")
    return
end

allowable = '[a-zA-Z][A-Za-z0-9_]*';
new_name = join(regexp(name,allowable, 'match'), '');

if ~isequal(new_name, name) 
    error(message('quantum:quantumCircuit:invalidName'));
end

end

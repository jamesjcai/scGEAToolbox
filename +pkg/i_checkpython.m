function [y] = i_checkpython

pe = pyenv;
y = ~isempty(pe.Version);

end
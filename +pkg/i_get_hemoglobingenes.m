function [g] = i_get_hemoglobingenes

pw1 = fileparts(mfilename('fullpath'));
txtfile = fullfile(pw1, '..','resources', 'HGNC', 'hemoglobin.txt');
warning off
t = readtable(txtfile, 'ReadVariableNames',false);
warning on
g = string(t.Var1);
% delete(fname);
end

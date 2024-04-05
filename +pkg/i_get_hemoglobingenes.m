function [g] = i_get_hemoglobingenes
%Get ribosomal genes

pw1 = fileparts(mfilename('fullpath'));
txtfile = fullfile(pw1, '..','resources', 'hemoglobin.txt');
warning off
t = readtable(txtfile, 'ReadVariableNames',false);
warning on
g = string(t.Var1);
% delete(fname);
end

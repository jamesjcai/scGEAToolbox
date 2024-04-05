function [g] = i_get_hemoglobingene
%Get ribosomal genes

pw1 = fileparts(mfilename('fullpath'));
txtfile = fullfile(pw1, '..','resources', ' hemoglobin.txt');
warning off
t = readtable(txtfile, 'Range', 'A:B');
warning on
g = string(t.ApprovedSymbol);
% delete(fname);
end

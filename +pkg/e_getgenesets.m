pw1 = fileparts(mfilename('fullpath'));
fname = fullfile(pw1, '..','resources', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
load(fname, 'T');
Ttfgn = T(T.mor > 0, :);
[gid, gnlist] = grp2idx(Ttfgn.target);
[tid, tflist] = grp2idx(Ttfgn.tf);
t = zeros(max(tid), max(gid));
t(sub2ind([max(tid), max(gid)], tid, gid)) = 1; % Ttfgn.mor;

setnames=string(tflist);
setgenes=string(gnlist);
setmatrx=logical(t);

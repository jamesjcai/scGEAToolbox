function [setmatrx, setnames, setgenes] = e_getgenesets(option)

if nargin<1, option = 'TF'; end

switch option
    case {'TF',1}

        pw1 = fileparts(mfilename('fullpath'));
        fname = fullfile(pw1, '..','resources', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
        load(fname, 'T');
        Ttfgn = T(T.mor > 0, :);
        [gid, gnlist] = grp2idx(Ttfgn.target);
        [tid, tflist] = grp2idx(Ttfgn.tf);
        t = zeros(max(tid), max(gid));
        t(sub2ind([max(tid), max(gid)], tid, gid)) = Ttfgn.mor;
        
        setmatrx=logical(t);
        setnames=string(tflist);
        setgenes=string(gnlist);
    case {2,'GSEA'}
        [~, ~, Col] = gui.i_selectMSigDBGeneSet('human',true);
        setnames = fields(Col);
        glist=[];
        for k=1:length(setnames)  
            glist=[glist;Col.(setnames{k}).geneSymbols];
        end
end



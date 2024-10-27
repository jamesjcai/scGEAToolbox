function [r] = sc_dissratio(X, genelist, vebrose)

    if nargin < 3, vebrose = true; end
    
    r = nan(width(X), 1);
    pw1 = fileparts(mfilename('fullpath'));
    dbfile1 = fullfile(pw1, '..', 'resources', 'scCancer', 'diss_genes.txt');
    if ~exist(dbfile1, 'file'), error('Missing file diss_genes.txt.'); end
    T = readtable(dbfile1, 'FileType', 'text', 'ReadVariableNames', false);
    dissog = string(T.Var1);
    
    idx = ismember(upper(genelist), upper(dissog));
    
    if sum(idx) > 0
        if vebrose
            fprintf('%d dissociation-associated genes found.\n', sum(idx));
        end
        lbsz = sum(X, 1);
        lbsz_mt = sum(X(idx, :), 1);
        r = transpose(lbsz_mt./lbsz);
    else
        if vebrose
            fprintf('No dissociation-associated genes found.\n');
        end
    end
end
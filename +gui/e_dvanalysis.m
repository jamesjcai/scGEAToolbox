function [T, X1, X2, g, xyz1, xyz2,...
    px1, py1, pz1, px2, py2, pz2] = e_dvanalysis(sce1, sce2, cL1, cL2)

if nargin<3
    cL1 = {'1'};
    cL2 = {'2'};
end
    if sce1.NumCells < 50 || sce2.NumCells < 50
        warning('One of groups contains too few cells (n < 50). The result may not be reliable.');
    end
    if sce1.NumGenes < 50 || sce2.NumGenes < 50
        warning('One of groups contains too few genes (n < 50). The result may not be reliable.');
    end

    if ~isequal(sce1.g, sce2.g)
        [g_ori, ia, ib] = intersect(sce1.g, sce2.g,'stable');
        X1_ori = sce1.X(ia, :);
        X2_ori = sce2.X(ib, :);
    else
        g_ori = sce1.g;
        X1_ori = sce1.X;
        X2_ori = sce2.X;
    end
    
    X1_ori = sc_norm(X1_ori,'type','libsize');
    X2_ori = sc_norm(X2_ori,'type','libsize');

    [T1, X1, g1, xyz1] = sc_splinefit(X1_ori, g_ori, true, false);
    [T1, idx1] = sortrows(T1,'genes','ascend');
    X1 = X1(idx1, :);
    g1 = g1(idx1);


    [T2, X2, g2, xyz2] = sc_splinefit(X2_ori, g_ori, true, false);
    [T2, idx2] = sortrows(T2,'genes','ascend');
    X2 = X2(idx2, :);
    g2 = g2(idx2);

    assert(isequal(g1, g2))
    g = g1;


    px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
    px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;

    %assignin("base","V1",[px1 py1 pz1]);
    %assignin("base","T1",T1);
    %assignin("base","xyz1",xyz1);

    v1=([px1 py1 pz1] - xyz1(T1.nearidx,:));
    v2=([px2 py2 pz2] - xyz2(T2.nearidx,:));

    DiffDist = vecnorm(v1 - v2, 2, 2);
    DiffSign = sign(vecnorm(v1,2,2)-vecnorm(v2,2,2));

    % Standardize distances
    ddz = zscore(DiffDist);
    % Assume Gaussian normal distribution for errors for p-values
    pval = 1 - normcdf(ddz);
    clear ddz;

    T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', cL1{1}));
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', cL2{1}));
    idxx = T1.(8)==1 | T2.(8)==1 | T1.(8) == max(T1.(8)) | T2.(8) == max(T2.(8));
    T1 = T1(:,1:end-1);
    T2 = T2(:,2:end);
    T2 = T2(:,1:end-1);

    T1.Properties.VariableNames{1} = 'gene';

    T = [T1 T2 table(DiffDist) table(DiffSign) table(pval)];

    T.DiffDist(idxx) = 0;
    T = sortrows(T,"DiffDist","descend");
    

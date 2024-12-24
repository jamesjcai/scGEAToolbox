function [T, X1, X2, g, xyz1, xyz2,...
    px1, py1, pz1, px2, py2, pz2] = e_dvanalysis_brennecke(sce1, sce2, cL1, cL2)

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
    
    [T1] = sc_hvg(X1_ori, g_ori, true, false);
    T1(T1.removedidx1 | T1.removedidx2,:)=[];
    [T2] = sc_hvg(X2_ori, g_ori, true, false);
    T2(T2.removedidx1 | T2.removedidx2,:)=[];
    [gene, idx1, idx2]=intersect(T1.genes, T2.genes);
    
    DiffDist = T1.residualcv2(idx1)-T2.residualcv2(idx2);
    DiffDistAbs = abs(DiffDist);
    DiffSign = sign(DiffDist);
    % Standardize distances
    ddz = zscore(DiffDistAbs);
    % Assume Gaussian normal distribution for errors for p-values
    pval = 1 - normcdf(ddz);
    % T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s', cL1{1}));
    % T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s', cL2{1}));
    % idxx = T1.(8)==1 | T2.(8)==1 | T1.(8) == max(T1.(8)) | T2.(8) == max(T2.(8));
    % T1 = T1(:,1:end-1);
    % T2 = T2(:,2:end);
    % T2 = T2(:,1:end-1);
    % 
    % T1.Properties.VariableNames{1} = 'gene';
    % 
    T = table(gene, DiffDist, DiffDistAbs, DiffSign, pval);
    T = sortrows(T,"DiffDistAbs","descend");
    

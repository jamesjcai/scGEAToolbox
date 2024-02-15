function [Tdv_final, Tig, influgenes] = sc_diffvar( X1, X2, genelist1, genelist2, fname, plotit )
    %{ 
    sc_variability computes the intersection of two experimets/batches and
    computes a cubic spline for each set to compute the variability.
    INPUT
        X ---------> Count matrix of set1 (Normalized)
        X2 --------> Count matrix of set2 (Normalized)
        genelist --> Gene list of set1
        genelist2 -> Gene list of set2
        fname -----> file name of tables
    OUTPUT
        Tdv_final -> Differentially variable table containing gene info
        Tig -------> Influenciable gene table containing polyfit info
        influgenes-> Influenciable gene list from intersected datasets
    NOTE
        1.- The intersection needs to be computed to use
            Tdv_final.dvgenes_idx as follows:
            [gene,~,~] = intersect(genelist, genelist2, 'stable'))
            gene(Tdv_final.dvgenes_idx)

        2.- Gene lists and tablers are sorted by genes to have same 
            indices in same places.

        3.- List is reduced by the 10% largest absolute differences 
            in polyfit-gene distance (dvgenes_idxf).
    USAGE
        Selected genes expressed in at least %7.5
        [X,genelist]=sc_selectg(X,genelist,1,0.075);
        [X]=sc_norm(X,'type','libsize');
        [X2,genelist2]=sc_selectg(X2,genelist2,1,0.075);
        [X2]=sc_norm(X2,'type','libsize');

        [Tdv, Tig, influgenes] = sc_variability( X, X2, genelist, genelist2, false);
    %}

    if nargin < 5, plotit = false; end

    % This needs to be done, otherwise polyfit gets junk
    sortit = true;

    [genelist3, irows, jrows] = intersect(genelist1, genelist2, 'stable');

    % T = table(lgu, lgcv, dropr, d, pval, fdr);
    % Compute cubic spline polynomial of noisy data
    [T1,~,genelist1,xyzp1] = sc_splinefit( X1(irows, :), genelist3, sortit);
    xyz1 = [T1.lgu, T1.lgcv, T1.dropr];
 
    % Compute cubic spline polynomial of noisy data
    [T2,~,genelist2,xyzp2]= sc_splinefit( X2(jrows, :), genelist3, sortit);
    xyz2 = [T2.lgu, T2.lgcv, T2.dropr];
    clear genelist3  irows jrows T1 T2;

    if plotit 
        figure;
        scatter3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), 'filled', 'MarkerFaceAlpha', .1);
        hold on
        plot3(xyzp1(:, 1), xyzp1(:, 2), xyzp1(:, 3), '-', 'linewidth', 4);
        hold on
        scatter3(xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), 'filled', 'MarkerFaceAlpha', .1);
        hold on
        plot3(xyzp2(:, 1), xyzp2(:, 2), xyzp2(:, 3), '-', 'linewidth', 4);
        xlabel('Mean, log');
        ylabel('CV, log');
        zlabel('Dropout rate (% of zeros)');
    end

    % Position vector referenced at polynomial
    % dsearchn(P,QP) utilizes P data points and QP query-test points (P>=QP)
    % xyz(idx) will sort the points to match xyzp order and size
    %[idx,dist] = dsearchn(xyz, xyzp); 
    [nearestidx1,dist1] = dsearchn(xyzp1, xyz1);

    Tdv = table(genelist1, dist1);
    Tdv.Properties.VariableNames = {'genes', 'dist'};
    [Tdv, idx] = sortrows(Tdv,'genes','descend');
    genelist1 = genelist1(idx);
    xyz1 = xyz1(idx,:);
    xyzp1 = xyzp1(idx,:); % Probably not needed
    clear dist1;

    % Position vector referenced at polynomial
    [nearestidx2,dist2] = dsearchn(xyzp2, xyz2);

    Tdv2 = table(genelist2, dist2);
    Tdv2.Properties.VariableNames = {'genes', 'dist'};
    [Tdv2, idx] = sortrows(Tdv2,'genes','descend');
    %genelist2 = genelist2(idx);
    %xyz2 = xyz2(idx,:);
    xyzp2 = xyzp2(idx,:); % Probably not needed
    clear dist2 genelist2 xyz2;

    % Take the 10% largest absolute differences in polyfit-gene distance
    abs_dist_diff = abs(Tdv.dist - Tdv2.dist);
    maxidx = round( length( abs_dist_diff )*0.1 );
    [abs_max_distf, dvgenes_idxf] = maxk( abs_dist_diff, maxidx);

    % Table containing distance metrics of differentially variable genes
    Tdv_final = table(Tdv.genes(dvgenes_idxf), Tdv.dist(dvgenes_idxf), Tdv2.dist(dvgenes_idxf), ...
                      abs_max_distf, dvgenes_idxf );

    clear Tdv Tdv2 idx;

    Tdv_final.Properties.VariableNames = {'genes', 'dist_set1','dist_set2','abs_diff','dvgenes_idx'};
   
    % Writting differential variability table 
    filenameT = '_diff_var.csv';
    filenameT = strcat(fname,filenameT);
    fprintf("File name is : %s",filenameT);
    writetable(Tdv_final,filenameT,'Delimiter',',');
    type(filenameT);

    % Termporarily shuted down
    active = false;
    if active 
    % Is polynomial difference saying something else? is this a scam?
    % kpolyg_idx contains the order for accessing xyzp (mapping)
    [~,distp] = dsearchn(xyzp1, xyzp2);

    % Look at top 10% farthest distances in polynomial xyzp order
    maxidx = round( length( distp )*.1);
    [max_dist, idx] = maxk( distp, maxidx);

    % Double mapping for obtaining genes in polynomial comparison 
    % Mapping data points to largest polynomial distance for obtain gene indices 
    npolyg_idx = dsearchn(xyz1(:,1:3), xyzp1(idx,1:3));

    % Genes for largest polynomial distances
    polygenes = genelist1(npolyg_idx);
    %polygenes = genelist(dvgenes_idx);

    % Table containing influenciable genes?
    Tig = table(polygenes, max_dist);
    Tig.Properties.VariableNames = {'genes', 'polyfit_dist'};

    % Writting differential variability polynomial fit table 
    %filenameT = 'diff_var_polyfit.csv';
    %writetable(Tig,filenameT,'Delimiter',',');
    %type 'diff_var_polyfit.csv';

    % Writting intersection of polyfit and differential variability
    [influgenes, ~, ~] = intersect(Tig.genes, Tdv_final.genes, 'stable');
    %filename = 'infugenes.csv';
    %writematrix(influgenes, filename);
    %type 'infugenes.csv';
    else 
        Tig = 0;
        influgenes = 0;
    end
end
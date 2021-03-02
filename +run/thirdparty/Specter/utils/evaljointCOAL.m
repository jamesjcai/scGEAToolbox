function results = evaljointCOAL(clusters,n1, n2, fea, fea_adt, n_clusters)
% clusters: each row is a clustering.     
    % clusters1 = clusters(1:n1, :);
    % [N, m] = size(clusters1);
    % X1 = zeros(m);
    % parfor i=1:1:m
    %     X1(i,:) = sum(clusters1 == clusters1(:,i));
    % end

    % clusters2 = clusters(n1+1:n1+n2, :);
    % [N, m] = size(clusters2);
    % X2 = zeros(m);
    % parfor i=1:1:m
    %     X2(i,:) = sum(clusters2 == clusters2(:,i));
    % end

    % %X = max(X1,X2);
    % if (n1 > 0)
    %     X1 = X1/max(max(X1));
    % end
    % if n2 > 0 
    %     X2 = X2/max(max(X2));
    % end
    % X = X1 + X2 + corrcoef(fea') + corrcoef(fea_adt');

    [N, m] = size(clusters);
    X = zeros(m);
    parfor i=1:1:m
        X(i,:) = sum(clusters == clusters(:,i));
    end
    X = 4*X/max(max(X)) + corrcoef(fea') + corrcoef(fea_adt');
    results = runLWEA(X, n_clusters);
end


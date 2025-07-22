% ensemble clustering using maximum AIR cluster
function results = evalEnsembleIdx(clusters, score)
% clusters: each row is a clustering.
[N, ~] = size(clusters);
X = zeros(N);
if strcmp(score, 'mi')
    for i = 1:1:N
        for j = 1:1:N
            X(i, j) = evalmutual(clusters(i, :), clusters(j, :));
        end
    end
elseif strcmp(score, 'nmi')
    for i = 1:1:N
        for j = 1:1:N
            X(i, j) = eval_nmi(clusters(i, :)', clusters(j, :)'); %eval_nmi, evalmutual,
        end
    end
elseif strcmp(score, 'ari')
    for i = 1:1:N
        for j = 1:1:N
            X(i, j) = eval_rand(clusters(i, :), clusters(j, :)); %eval_nmi, evalmutual,
        end
    end
elseif strcmp(score, 'ri')
    for i = 1:1:N
        for j = 1:1:N
            X(i, j) = rand_index(clusters(i, :), clusters(j, :)); %eval_nmi, evalmutual,
        end
    end
else
    disp('score is undefined');
end

sumX = sum(X);
[~, idx] = max(sumX);
results = clusters(idx, :)';
end

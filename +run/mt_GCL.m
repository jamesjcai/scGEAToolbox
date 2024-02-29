function [v] = mt_GCL(X, k)

if nargin < 2, k = 50; end

[v] = gcl_ori(X, k);

end


function [gcl_output] = gcl_ori(data, num_divisions)
% Calculate the GCL of the input data with bootstrap
% Guy Amit, guy1.amit@gmail.com, Orr Levy, Dana Vaknin, Tom Snir, Sol
% Efroni, Peter Castaldi, Yang-Yu Liu, Haim Cohen, Amir Bashan.
% Based on Bias Corrected Distance Correlation Szekely, G. J., & Rizzo,
% M. L. (2013). The distance correlation t-test of independence in
% high dimension. Journal of Multivariate Analysis, 117, 193-213.
%
% data - Input data such that [num_genes, num_cells] = size(data)
% num_division - Number of random gene division for calculation.
% gcl_output - The GCL of the data

% https://github.com/guy531/gcl/blob/master/gcl.m

[num_genes, n] = size(data);
gcl_output = zeros(num_divisions, 1);

% Creating functions
Vn = @(Aij, Bij) (1 / (n * (n - 3))) .* (sum(sum(Aij.*Bij)) - (n / (n - 2)) * diag(Aij)' * diag(Bij));
Rn = @(Aij, Bij) Vn(Aij, Bij) ./ sqrt(Vn(Aij, Aij).*Vn(Bij, Bij));

rng("shuffle");
for i = 1:num_divisions
    % Random divisions
    random_genes = randperm(num_genes);
    X1 = data(random_genes(1:floor(num_genes/2)), :)';
    X2 = data(random_genes(floor(num_genes/2)+1:end), :)';

    % Aij matrices
    d1 = pdist2(X1, X1);
    m1 = mean(d1);
    M1 = mean(d1(:));
    Aij1 = d1 - m1' * ones(1, n);
    Aij1 = Aij1 - ones(n, 1) * m1;
    Aij1 = Aij1 + M1;
    Aij1 = Aij1 - d1 / n;
    Aij1(1:n+1:end) = m1 - M1;
    Aij1 = (n / (n - 1)) * Aij1;


    d2 = pdist2(X2, X2);
    m2 = mean(d2);
    M2 = mean(d2(:));
    Aij2 = d2 - m2' * ones(1, n);
    Aij2 = Aij2 - ones(n, 1) * m2;
    Aij2 = Aij2 + M2;
    Aij2 = Aij2 - d2 / n;
    Aij2(1:n+1:end) = m2 - M2;
    Aij2 = (n / (n - 1)) * Aij2;

    % Calculating bcdcorr
    gcl_output(i) = Rn(Aij2, Aij1);

end


end
function [A] = sc_pcnet(X, ncom, fastersvd, dozscore, guiwaitbar)
%Construct GRN network A using PC regression (pcnet)
%
% [X]=log(1+sc_norm(X));     % pcnet input should be LogNormalized
% [A]=sc_pcnet(X,ncom);      % X = expression matrix of genes x cells
% ncom - number of components used (default=3)
% ref: https://rdrr.io/cran/dna/man/PCnet.html
% https://github.com/cran/dna/blob/master/src/rpcnet.c
% https://rdrr.io/cran/dna/f/inst/doc/Introduction.pdf

arguments
    X double
    ncom(1, 1) {mustBeNumeric} = 3
    fastersvd(1, 1) logical = false
    dozscore(1, 1) logical = false
    guiwaitbar(1, 1) logical = false
end

%if nargin<5 || isempty(guiwaitbar), guiwaitbar=false; end
%if nargin<4 || isempty(dozscore), dozscore=true; end
%if nargin<3 || isempty(fastersvd), fastersvd=false; end
%if nargin<2 || isempty(ncom), ncom=3; end


opts.maxit = 150;

% if fastersvd
%     pw1 = fileparts(mfilename('fullpath'));
%     pth = fullfile(pw1, '+run', 'thirdparty', 'faster_svd', 'lmsvd');
%     if ~(ismcc || isdeployed), addpath(pth); end
% end

% LogNormalize: Feature counts for each cell are divided by
% the total counts for that cell and multiplied by the
% scale.factor. This is then natural-log transformed using log1p.
% https://satijalab.org/seurat/reference/normalizedata
%[X]=sc_norm(X);

X = X.';
if dozscore
    % whos("X")
    X = zscore(X);
end
n = size(X, 2);
A = 1 - eye(n);

if guiwaitbar
    fw = gui.gui_waitbar_adv;
end
for k = 1:n
    % fprintf('...... %d/%d\n',k,n);
    if guiwaitbar
        gui.gui_waitbar_adv(fw, k/n);
    end
    y = X(:, k);
    Xi = X;
    Xi(:, k) = [];
    if fastersvd
        [~, ~, coeff] = lmsvd(Xi, ncom, opts);
    else
        [~, ~, coeff] = svds(Xi, ncom);
        % https://www.mathworks.com/matlabcentral/fileexchange/47835-randomized-singular-value-decomposition
        %[~,~,coeff]=rsvd(Xi,ncom);   % not recommanded
    end
    score = Xi * coeff;
    % [coeff,score]=pca(Xi);
    % coeff=coeff(:,1:ncom);
    score = (score ./ (vecnorm(score).^2));
    Beta = sum(y.*score);
    A(k, A(k, :) == 1) = coeff * Beta';
end
if guiwaitbar, gui.gui_waitbar_adv(fw); end
end

% [PCALoadings,PCAScores] = pca(Xi,"NumComponents",ncom);
% betaPCR = regress(y-mean(y), PCAScores(:,1:2));
% https://www.mathworks.com/help/stats/partial-least-squares-regression-and-principal-components-regression.html?prodcode=ST&language=en
% https://blogs.sas.com/content/iml/2017/10/25/principal-component-regression-drawbacks.html

%{
library(dna)

X1=rbind(
c(2.5,6.7,4.5,2.3,8.4,3.1),
c(1.2,0.7,4.0,9.1,6.6,7.1),
c(4.3,-1.2,7.5,3.8,1.0,9.3),
c(9.5,7.6,5.4,2.3,1.1,0.2))
s=PCnet(X1)
print(round(s,4))

# small example using PCnet with 2 principal components,
# data rescaled, and scores symmetrized and rescaled
s2=PCnet(X1,ncom=3,rescale.data=TRUE,symmetrize.scores=FALSE,rescale.scores=FALSE)
print(round(s2,4))
library(dna)
X=rbind(
c(2.5,6.7,4.5,2.3,8.4,3.1),
c(1.2,0.7,4.0,9.1,6.6,7.1),
c(4.3,-1.2,7.5,3.8,1.0,9.3),
c(9.5,7.6,5.4,2.3,1.1,0.2))
s=PCnet(X,ncom=3,rescale.data=TRUE,symmetrize.scores=FALSE)
print(round(s,4))

%}

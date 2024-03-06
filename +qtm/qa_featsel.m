%{
load GSM3308547_GSM330854.mat
y = sce.list_cell_attributes{2};
X = log(1+sc_norm(sce.X));

K = 10;
nbins=20;
probs = [0,(1:(nbins-1))/nbins,1];
edges = quantile(y,probs);
y = discretize(y,edges,IncludedEdge="right");
X = 1+(X>0);

[p,n] = size(X);
nbinsX=repmat(2,p,1);
R0 = iComputeMIXX(X',nbinsX);



J = iComputeMIXY(X',y,nbinsX,nbinsY);
R = R0/(K-1);
%}
K=5
fun = @(alpha)howmany(alpha,R,J) - K;
alphasol = fzero(fun,1);

[~,xsol] = howmany(alphasol,R,J);


sce.g(find(xsol.BestX))



% rng default
% mdl = fitrtree(X,Y,CrossVal="on",Holdout=0.2);
% kfoldLoss(mdl)
% 
% rng default
% X2 = X(:,find(xsol.BestX));
% mdl2 = fitrtree(X2,Y,CrossVal="on",Holdout=0.2);
% kfoldLoss(mdl2)


% =================================


% function [N,p,X,y] = iSyntheticData1
% rng default
% N = 10000;
% p = 30;
% useful = [6,8,9,11,13];
% C = randn(p,p);
% R = corrcov(C'*C);
% X = mvnrnd(zeros(p,1),R,N);
% % Make features 15 to 19 highly correlated with useful features:
% % 15 -> 6
% % 16 -> 8
% % 17 -> 9
% % 18 -> 11
% % 19 -> 13
% corrStd = 0.1;
% X(:,15:19) = X(:,useful) + corrStd*randn(N,5);
% noiseStd = 0.1;
% t = 0.5*cos(X(:,11)) + sin(X(:,9).*X(:,8)) + 0.5*X(:,13).*X(:,6) + noiseStd*randn(N,1);
% y = rescale(t,0,1);
% X = zscore(X);
% end

%This code creates the iBinXY helper function. Note that this helper function uses the iBinPredictors helper function.

% function [binnedXInfo,binnedYInfo] = iBinXY(X,Y,nbins)
% binnedXInfo = iBinPredictors(X,nbins);
% binnedYInfo = iBinPredictors(Y,nbins);
% end

%This code creates the iComputeMIXX helper function. Note that this helper function uses the iComputeMIXIXJ helper function.

function Q = iComputeMIXX(binnedX,nbinsX)
p = size(binnedX,2);
Q = zeros(p,p);
for i = 1:p
    for j = i+1:p
        Xi = binnedX(:,i);
        Xj = binnedX(:,j);
        nbinsXi = nbinsX(i);
        nbinsXj = nbinsX(j);
        Q(i,j) = iComputeMIXIXJ(Xi,Xj,nbinsXi,nbinsXj);
    end
end
Q = Q + Q';
end

%This code creates the iComputeMIXIXJ helper function.

function mi = iComputeMIXIXJ(Xi,Xj,nbinsXi,nbinsXj)
N = size(Xi,1);
NXY = zeros(nbinsXi,nbinsXj);
for n = 1:N
    a = Xi(n);
    b = Xj(n);
    NXY(a,b) = NXY(a,b) + 1;
end
NX = sum(NXY,2);
NY = sum(NXY,1);
NXNYDivN = NX.*NY/N;
nonzeroID = NXY ~= 0;
nonzeroNXY = NXY(nonzeroID);
nonzeroNXNYDivN = NXNYDivN(nonzeroID);
mi = sum((nonzeroNXY./N).*log(nonzeroNXY./nonzeroNXNYDivN),'all');
end
%This code creates the iComputeMIXY helper function. Note that this helper function uses the iComputeMIXIXJ helper function.

function f = iComputeMIXY(binnedX,binnedY,nbinsX,nbinsY)
p = size(binnedX,2);
f = zeros(p,1);
for i = 1:p
    Xi = binnedX(:,i);
    nbinsXi = nbinsX(i);
    f(i) = iComputeMIXIXJ(Xi,binnedY,nbinsXi,nbinsY);
end
end
%This code creates the howmany helper function.

function [n,result] = howmany(alpha,R,J)
Q = qubo((1-alpha)*R - alpha*diag(J));
result = solve(Q);
n = sum(result.BestX);
end

%This code creates the iBinPredictors helper function. Note that this helper function uses the iDiscretize helper function.

function binnedInfo = iBinPredictors(X,nbins)
[N,p] = size(X);
binnedX = zeros(N,p);
edgesX = cell(1,p);
iscatX = false(1,p);
nbinsX = zeros(1,p);
istbl = istable(X);
for i = 1:p
    if istbl
        oneX = X{:,i};
    else
        oneX = X(:,i);
    end
    nbinsi = min(nbins, numel(unique(oneX)));
    [binnedX(:,i),edgesX{i},iscatX(i),nbinsX(i)] = iDiscretize(oneX,nbinsi);
end
binnedInfo = struct;
binnedInfo.binned = binnedX;
binnedInfo.edges = edgesX;
binnedInfo.iscat = iscatX;
binnedInfo.nbins = nbinsX;
end

%This code creates the iDiscretize helper function.

function [binnedx,edges,iscat,nbins] = iDiscretize(x,nbins)
if iscategorical(x)
    binnedx = double(x);
    edges = [];
    iscat = true;
    nbins = numel(unique(x));
else
    probs = [0,(1:(nbins-1))/nbins,1];
    edges = quantile(x,probs);
    binnedx = discretize(x,edges,IncludedEdge="right");
    iscat = false;
end
end

function [y, S, F, ydata,alphaK,timeOurs,converge,LF] = SIMLR_pearson(X, c, k, ifimpute,normalize)

%%%
if nargin==2
    k=10;
    ifimpute = 0;
    normalize = 0;
end

if nargin==3
    ifimpute = 0;
    normalize = 0;
end

if nargin==4
    normalize = 0;
end


if ifimpute
    X = X';
    [I,J] = find(X==0);
    Xmean = mean(X);
    X(sub2ind(size(X),I,J)) = Xmean(J);
    X = X';
end

if normalize
    X = X';
    X = X - min(X(:));
    X = X / max(X(:));
    X = bsxfun(@minus, X, mean(X, 1));
    X = X';
end

t0 = tic;
order = 2;
no_dim=c;


NITER = 30;
num = size(X,1);
r = -1;
beta = 0.8;
D_Kernels = multipleK_pearson(X);

% https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby076/5077112#122849978


clear X
alphaK = 1/(size(D_Kernels,3))*ones(1,size(D_Kernels,3));
distX = mean(D_Kernels,3);
[distX1, idx] = sort(distX,2);
A = zeros(num);
di = distX1(:,2:(k+2));
rr = 0.5*(k*di(:,k+1)-sum(di(:,1:k),2));
id = idx(:,2:k+2);
temp = (repmat(di(:,k+1),1,size(di,2))-di)./repmat((k*di(:,k+1)-sum(di(:,1:k),2)+eps),1,size(di,2));
a = repmat((1:num)',1,size(id,2));
A(sub2ind(size(A),a(:),id(:)))=temp(:);
if r <= 0
    r = mean(rr);
end
lambda = max((mean(rr)),0);
A(isnan(A))=0;
S0 = max(max(distX))-distX;
S0 = Network_Diffusion(S0,k);
S0 = NE_dn(S0,'ave');
S= (S0 + S0')/2;
D0 = diag(sum(S,order));
L0= D0-S;
[F, ~, evs]=eig1(L0, c, 0);
F = NE_dn(F,'ave');
for iter = 1:NITER
    distf = L2_distance_1(F',F');
    A = zeros(num);
    %b = idx(:,2:(2*k+2));
    b = idx(:,2:end);
    a = repmat((1:num)',1,size(b,2));
    inda = sub2ind(size(A),a(:),b(:));
    ad = reshape((distX(inda)+lambda*distf(inda))/2/r,num,size(b,2));
    ad = projsplx_c(-ad')';
    A(inda) = ad(:);
    A(isnan(A))=0;
    S = (1-beta)*A+beta*S;
    S = Network_Diffusion(S,k);
    S= (S + S')/2;
    D = diag(sum(S,order));
    L = D - S;
    F_old = F;
    [F, ~, ev]=eig1(L, c, 0);
    F = NE_dn(F,'ave');
    F = (1-beta)*F_old+beta*F;
    evs(:,iter+1) = ev;
    for i = 1:size(D_Kernels,3)
        temp = (eps+D_Kernels(:,:,i)).*(eps+S);
        DD(i) = mean(sum(temp));
    end
    alphaK0 = umkl_bo(DD);
    alphaK0 = alphaK0/sum(alphaK0);
    alphaK = (1-beta)*alphaK + beta*alphaK0;
    alphaK = alphaK/sum(alphaK);
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    converge(iter) = fn2-fn1;
    if iter<10
        if (ev(end) > 0.000001)
            lambda = 1.5*lambda;
            r = r/1.01;
        end
    else
        if (converge(iter)>1.01*converge(iter-1))
            S = S_old;
            if converge(iter-1) > 0.2
                warning('Maybe you should set a larger value of c');
            end
            break;
        end
    end
    S_old = S;
    distX = Kbeta(D_Kernels,alphaK');
    [~, idx] = sort(distX,2);
end
LF = F;
D = diag(sum(S,order));
L = D - S;
[U,~] = eig(L);
if length(no_dim)==1
    F = tsne_p_bo((S),[], U(:,1:no_dim));
else
    clear F;
    for i = 1:length(no_dim)
        F{i} = tsne_p_bo((S),[], U(:,1:no_dim(i)));
    end
end
timeOurs = toc(t0);

[~,center] = litekmeans(LF, c,'replicates',20);
[~,center] = min(dist2(center,LF),[],2);
y = litekmeans(F,c,'Start',center);
ydata = tsne_p_bo(S);
end
function thisP = umkl_bo(D,beta)
if nargin<2
    beta = 1/length(D);
end
tol = 1e-4;
u = 20;logU = log(u);
[H, thisP] = Hbeta(D, beta);
betamin = -Inf;
betamax = Inf;
% Evaluate whether the perplexity is within tolerance
Hdiff = H - logU;
tries = 0;
while (abs(Hdiff) > tol) && (tries < 30)
    
    % If not, increase or decrease precision
    if Hdiff > 0
        betamin = beta;
        if isinf(betamax)
            beta = beta * 2;
        else
            beta = (beta + betamax) / 2;
        end
    else
        betamax = beta;
        if isinf(betamin)
            beta = beta / 2;
        else
            beta = (beta + betamin) / 2;
        end
    end
    
    % Recompute the values
    [H, thisP] = Hbeta(D, beta);
    Hdiff = H - logU;
    tries = tries + 1;
end
end



function D_Kernels = multipleK_pearson(x)


N = size(x,1);
KK = 0;
sigma = [2:-0.25:1];
Diff = (dist2_pearson(x));
[T,~]=sort(Diff,2);
[~,n]=size(Diff);
allk = 10:2:30;
t=1;
for l = 1:length(allk)
    if allk(l) < (size(x,1)-1)
        TT=mean(T(:,2:(allk(l)+1)),2)+eps;
        Sig=(repmat(TT,1,n)+repmat(TT',n,1))/2;
        Sig=Sig.*(Sig>eps)+eps;
        for j = 1:length(sigma)
            W=normpdf(Diff,0,sigma(j)*Sig);
            Kernels(:,:,KK+t) = (W + W')/2;
            t = t+1;
        end
    end
end

for i = 1:size(Kernels,3)
    K = Kernels(:,:,i);
    k = 1./sqrt(diag(K)+1);
    %G = K.*(k*k');
    G = K;
    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
  
end

end
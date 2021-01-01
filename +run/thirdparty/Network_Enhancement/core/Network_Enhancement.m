function [ W_out] = Network_Enhancement( W_in,order,K,alpha)
%%%W_in the the input network of size N X N
%%%K the number of neighbors
%%% alpha is the regularization parameter
%%% order determines the extent of diffusion. Typical values are 0.5,1,2.



if nargin<2
    K =min(20,ceil(length(W_in)/10));
    alpha = 0.9;
    order = 2;
end

if nargin<3
    K =min(20,ceil(length(W_in)/10));
    alpha = 0.9;
end

if nargin<4
    alpha = 0.9;
end

W_in1 = W_in.*(1-eye(length(W_in)));
zeroindex = find(sum(abs(W_in1))>0);
W0 = W_in(zeroindex,zeroindex);
W = NE_dn(W0,'ave');
W = (W+W')/2;

DD = sum(abs(W0));

if length(unique(W(:)))==2
    P = W;
else
    P = (dominateset(double(abs(W)),min(K,length(W)-1))).*sign(W);
end
P = P + (eye(length(P))+diag(sum(abs(P'))));
%P = P + (eye(length(P)));
P = TransitionFields(P);
[U,D] = eig(full(P));
d = real(diag(D-eps));
d = (1-alpha)*d./(1-alpha*d.^order);
D = diag(real(d));
W = U*D*U';
W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),(DD));
W=D*(W);
W(W<0)=0;
W = (W+W')/2;
W_out = zeros(size(W_in));
W_out(zeroindex,zeroindex) = W;
end


function [W,H] = nndsvd(A,k,flag)
%
% This function implements the NNDSVD algorithm described in [1] for
% initializattion of Nonnegative Matrix Factorization Algorithms.
%
% [W,H] = nndsvd(A,k,flag);
%
% INPUT
% ------------
%
% A    : the input nonnegative m x n matrix A
% k    : the rank of the computed factors W,H
% flag : indicates the variant of the NNDSVD Algorithm
%        flag=0 --> NNDSVD
%        flag=1 --> NNDSVDa
%        flag=2 --> NNDSVDar
%
% OUTPUT
% -------------
%   
% W   : nonnegative m x k matrix
% H   : nonnegative k x n matrix
%
% 
% References:
% 
% [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
%     start for nonnegative matrix factorization, Pattern Recognition,
%     Elsevier
%
% This code is kindly provided by the authors for research porpuses.
% - Efstratios Gallopoulos (stratis@ceid.upatras.gr)
% - Christos Boutsidis (boutsc@cs.rpi.edu)     
%

%check the input matrix
if numel(find(A<0)) > 0
    error('The input matrix contains negative elements !')
end

%size of input matrix
[m,n] = size(A);

%the matrices of the factorization
W = zeros(m,k);
H = zeros(k,n);

%1st SVD --> partial SVD rank-k to the input matrix A. 
[U,S,V] = svds(A,k);

%-------------------------------------------------------
% We also recommend the use of propack for the SVD
% 1st SVD --> partial SVD rank-k ( propack )
% OPTIONS.tol  = 0.00001;               % remove comment to this line
% [U,S,X] = LANSVD(A,k,'L',OPTIONS);    % remove comment to this line 
%-------------------------------------------------------

%choose the first singular triplet to be nonnegative
W(:,1)     =  sqrt(S(1,1)) * abs(U(:,1) );         
H(1,:)     =  sqrt(S(1,1)) * abs(V(:,1)'); 

% second SVD for the other factors (see table 1 in our paper)
for i=2:k
    uu = U(:,i); vv = V(:,i);
    uup = pos(uu); uun = neg(uu) ;
    vvp = pos(vv); vvn = neg(vv);
    n_uup = norm(uup);
    n_vvp = norm(vvp) ;
    n_uun = norm(uun) ;
    n_vvn = norm(vvn) ;
    termp = n_uup*n_vvp; termn = n_uun*n_vvn;
    if (termp >= termn)
        W(:,i) = sqrt(S(i,i)*termp)*uup/n_uup; 
        H(i,:) = sqrt(S(i,i)*termp)*vvp'/n_vvp;
    else
        W(:,i) = sqrt(S(i,i)*termn)*uun/n_uun; 
        H(i,:) = sqrt(S(i,i)*termn)*vvn'/n_vvn;
    end
end
%------------------------------------------------------------

%actually these numbers are zeros
W(find(W<0.0000000001))=0;
H(find(H<0.0000000001))=0;

%fill in the zero elements with the average : NNDSVDa
if flag==1
   ind1      =  find(W==0) ;
   ind2      =  find(H==0) ;
   average   =  mean(A(:)) ; 
   W( ind1 ) =  average    ; 
   H( ind2 ) =  average    ;

% fill in the zero elements with random values in the space :[0:average/100]
% NNDSVDar
elseif flag==2
   ind1      =  find(W==0) ;
   ind2      =  find(H==0) ;
   n1        =  numel(ind1);
   n2        =  numel(ind2);
   
   average   =  mean(A(:))       ;
   W( ind1 ) =  (average*rand(n1,1)./100)  ; 
   H( ind2 ) =  (average*rand(n2,1)./100)  ;   
end
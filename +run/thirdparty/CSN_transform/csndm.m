function ndm = csndm(data,alpha,boxsize,normalize)
%Construction of network degree matrix
%The function performs the transformation from gene expression matrix to
%network degree matrix (ndm).
%data: Gene expression matrix (TPM/RPKM/FPKM/count), rows = genes, columns = cells
%alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01
%boxsize: Size of neighborhood, Default = 0.1 (nx(k) = ny(k) = 0.1*n)
%normalize: 1  result is normalized (Default); 0  result is not normalized           
 
if nargin < 4 || isempty(normalize)
    normalize = 1;
end
if nargin < 3 || isempty(boxsize)
    boxsize = 0.1;
end
if nargin <2 || isempty(alpha)
    alpha = 0.01;
end
 
%Define the neighborhood of each plot
[n1,n2] = size(data);
upper = zeros(n1,n2);
lower = zeros(n1,n2);
for i = 1 : n1
    [s1,s2] = sort(data(i,:));
    n0 = n2-sum(sign(s1));
    h = round(boxsize/2*sum(sign(s1)));
    k = 1;
    while k <= n2
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)
            s = s+1;
        end
        if s >= h
            upper(i,s2(k:k+s)) = data(i,s2(k));
            lower(i,s2(k:k+s)) = data(i,s2(k));
        else
            upper(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
            lower(i,s2(k:k+s)) = data(i,s2(max(n0*(n0>h)+1,k-h)));
        end
        k = k+s+1;
    end
end
 
%If gene expression matrix is sparse, use the sparse matrix will accelerate
%the calculation and reduce memory footprint
%data = sparse(data); upper = sparse(upper); lower = sparse(lower);
 
%Construction of network degree matrix
ndm = zeros(n1,n2);
B = zeros(n1,n2);
p = -icdf('norm',alpha,0,1);
for k = 1 : n2
    for j = 1 : n2
        B(:,j) = data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k) & data(:,k);
    end
    %B = sparse(B);
    a = sum(B,2);
    csn = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
    csn = (csn > p);
    ndm(:,k) = sum(csn,2) - diag(csn);
    disp(['Cell ' num2str(k) ' is completed']);
end
 
%Normalization of network degree matrix
if normalize
    ndm = bsxfun(@rdivide,ndm,sum(ndm))*mean(sum(sign(ndm)))^2/2000;
end

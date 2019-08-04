function csn = csnet(data,c,alpha,boxsize,weighted)
%Construction of cell-specific networks
%The function performs the transformation from gene expression matrix to
%cell-specific network (csn).
%data: Gene expression matrix, rows = genes, columns = cells
%c: Construct the CSNs for all cells, set c = [] (Default);
%   Construct the CSN for cell k, set  c = k
%alpha: Significant level (eg. 0.001, 0.01, 0.05 ...)
%       larger alpha leads to more edges, Default = 0.01
%boxsize: Size of neighborhood, Default = 0.1
%weighted: 1  edge is weighted
%          0  edge is not weighted (Default)
%csn: Cell-specific network, the kth CSN is in csn{k}
%     rows = genes, columns = genes

% Too many cells or genes may lead to out of memory.

 
if nargin < 5 || isempty(weighted)
    weighted = 0;
end
if nargin < 4 || isempty(boxsize)
    boxsize = 0.1;
end
if nargin <3 || isempty(alpha)
    alpha = 0.01;
end
 
[n1,n2] = size(data);
if nargin <2 || isempty(c)
    c = 1 : n2;
end
 
%Define the neighborhood of each plot
upper = zeros(n1,n2);
lower = zeros(n1,n2);
for i = 1 : n1
    [s1,s2] = sort(data(i,:));
    n3 = n2-sum(sign(s1));
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
            lower(i,s2(k:k+s)) = data(i,s2(max(n3*(n3>h)+1,k-h)));
        end
        k = k+s+1;
    end
end
 
%Construction of cell-specific network
csn = cell(1,n2);
B = zeros(n1,n2);
p = -icdf('norm',alpha,0,1);
for k = c
    for j = 1 : n2
        B(:,j) = data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k);
    end
    a = sum(B,2);
    d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
    d(1 : n1+1 : end) = 0;
    if weighted
        csn{k} = d.*(d > 0);
    else
        csn{k} = sparse(d > p);
    end
    disp(['Cell ' num2str(k) ' is completed']);
end

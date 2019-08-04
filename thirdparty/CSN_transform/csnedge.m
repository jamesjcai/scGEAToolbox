function edge = csnedge(gx,gy,boxsize)
%The normalized statistic of edge x-y
%gx gy: Gene expression values of gene x and gene y
%       If there are n cells, gx and gy are 1-by-n vectors
%boxsize: Size of neighborhood, Default = 0.1
%edge: 1-by-n vector, the normalized statistic of edge x-y in all cells

if nargin < 3
    boxsize = 0.1;
end
 
%Define the neighborhood of each plot
n = length(gx);
upper = zeros(1,n);
lower = zeros(1,n);
a = zeros(2,n);
B = cell(1,2);
for i = 1 : 2
    g = gx*(i==1)+gy*(i==2);
    [s1,s2] = sort(g);
    n0 = n-sum(sign(s1));
    h = round(boxsize/2*sum(sign(s1)));
    k = 1;
    while k <= n
        s = 0;
        while k+s+1 <= n && s1(k+s+1) == s1(k)
            s = s+1;
        end
        if s >= h
            upper(s2(k:k+s)) = g(s2(k));
            lower(s2(k:k+s)) = g(s2(k));
        else
            upper(s2(k:k+s)) = g(s2(min(n,k+s+h)));
            lower(s2(k:k+s)) = g(s2(max(n0*(n0>h)+1,k-h)));
        end
        k = k+s+1;
    end
    
    B{i} = bsxfun(@le,g',upper) & bsxfun(@ge,g',lower);
    a(i,:) = sum(B{i});
end
 
%Calculate the normalized statistic of edge x-y in all cells
edge = (sum(B{1} & B{2})*n-a(1,:).*a(2,:))./sqrt(a(1,:).*a(2,:) ...
    .*(n-a(1,:)).*(n-a(2,:))/(n-1));

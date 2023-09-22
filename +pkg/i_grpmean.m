function [a, fano, cv] = i_grpmean(x, c)

n = numel(unique(c));
a = zeros(size(x, 1), n);
%atrim=zeros(size(x,1),n);
fano = zeros(size(x, 1), n);
cv = zeros(size(x, 1), n);
for k = 1:n
    %    atrim(:,k)=full(trimmean(x(:,c==k),10,2));
    a(:, k) = full(mean(x(:, c == k), 2));
    fano(:, k) = full(var(x(:, c == k), 0, 2));
    cv(:, k) = full(std(x(:, c == k), 0, 2));
end
if nargout > 1
    fano = fano ./ a;
end
if nargout > 2
    cv = cv ./ a;
end
%if nargout==1
%    a=atrim;
%end
end

function x = tridiag(a,b,c,d) 
% Function tridiag:
%    Inverts a tridiagonal system whose lower, main and upper diagonals
%    are respectively given by the vectors a, b and c. d is the right­hand
%    side. The result is placed in x.


[N,dim] = size(d);
assert(length(a) == N-1 && length(c) == N-1);
assert(length(b) == N);
x = zeros(N,dim);

c_new = zeros(N-1,1);
d_new = zeros(N,dim);

% Eliminate the lower diagonal

c_new(1) = c(1)/b(1);
d_new(1,:) = d(1,:)/b(1);

for i=2:N-1
    c_new(i) = c(i)/(b(i)-a(i-1)*c_new(i-1));
    d_new(i,:) = (d(i,:)-a(i-1)*d_new(i-1,:))/(b(i)-a(i-1)*c_new(i-1));  
end

d_new(N,:) = (d(N,:)-a(N-1)*d_new(N-1,:))/(b(N)-a(N-1)*c_new(N-1));

% Perform the backsubstitution

x(N,:) = d_new(N,:);

for j = 1:(N-1)
    k = N-j ;
    x(k,:) = d_new(k,:)-c_new(k)*x(k+1,:);
end

end



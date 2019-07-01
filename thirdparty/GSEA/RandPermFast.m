function p = RandPermFast(n)
% Random permutation using Fisher-Yates-Knuth Shuffle Algorithm

p = 1:n; % Identity permutation
for k = n:-1:2
    r = 1 + floor(rand*k); % random integer between 1 and k-m
    t = p(r);
    p(r) = p(k); % Swap(p(r),p(k))
    p(k) = t;
end
% Author: Modified package by Van Hoan Do
function D2 = mydist(x,y)
% compute the similarity distance between x and y
    % D2 = sum(x' == y');
    [m,p] = size(y);
    D2 = zeros(m,1);
    for q = 1:p
        D2 = D2 + (x(q) == y(:,q)); 
    end
end

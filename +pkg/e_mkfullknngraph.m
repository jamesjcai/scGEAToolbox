function [G, W] = e_mkfullknngraph(s, k)
if nargin < 2, k = 5; end

[A, W] = sc_knngraph(s, k);
G = graph(A);
b = conncomp(G);
n = max(b);
D = zeros(n, n);
Didx = cell(n, n);
for k = 1:n - 1
    for l = k:n
        i1 = find(b == k);
        i2 = find(b == l);
        d = pdist2(s(i1, :), s(i2, :));
        [D(k, l), idx] = min(d(:));
        [x, y] = ind2sub(size(d), idx);
        Didx{k, l} = [i1(x), i2(y)];
        %Didx{l,k}=Didx{k,l};
    end
end
DG = graph(D+D');
T = minspantree(DG);
Te = T.Edges;
for k = 1:size(Te, 1)
    x = Te.EndNodes(k, 1);
    y = Te.EndNodes(k, 2);
    a = Didx{x, y}(1);
    b = Didx{x, y}(2);
    A(a, b) = 1;
    W(a, b) = norm(s(a, :)-s(b, :));
    W(b, a) = W(a, b);
    G = addedge(G, a, b, 1.0);
end
%A=G.adjacency;

end

%{
figure;
hold on
for i = 1 : size(A,2)
    for j = 1 : size(A,1)
        if A(i,j)>0
            if size(s,2)>=3
                line([s(i,1),s(j,1)],...
                    [s(i,2),s(j,2)],...
                    [s(i,3),s(j,3)],...
                    'Color','red');
            else
                line([s(i,1),s(j,1)],...
                    [s(i,2),s(j,2)],'Color','red');
            end
        end
    end
end
hold off
%}
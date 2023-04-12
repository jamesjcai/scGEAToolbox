function [Y,idxNN] = isomap(X,K,d)
if nargin<3, d=3; end
if nargin<2, K=4; end

% https://github.com/gionuno/isomap
    Tree = createns(X,'NSMethod','kdtree');
    
    idxNN = knnsearch(Tree,X,'K',K+1);
    idxNN = idxNN(:,2:end);
    
    N = size(X,1);
    
    W = inf*ones(N,N);
    A = zeros(N,N);
    
    for i = 1:N
        %disp(i);
        for k = 1:K
            W(i,idxNN(i,k)) = norm(X(i,:)-X(idxNN(i,k),:));
            A(i,idxNN(i,k)) = 1;
            A(idxNN(i,k),i) = 1;
        end
    end
    
    D = floyd_warshall(min(W,W'));
    J = eye(N)-ones(N,N)/N;
    [Y,E] = eig(-0.5*J*(D.^2)*J);
    Y = Y(:,2:d+1);
end


function D = floyd_warshall(D_0)
    D = D_0;
    N = size(D,1);
    for k = 1:N
        disp(k);
        D = min(D,repmat(D(:,k),[1,N])+repmat(D(k,:),[N,1]));
    end
end
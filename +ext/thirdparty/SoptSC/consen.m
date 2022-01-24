function [all_eigs,M_all] = consen(S,idx,tau)

flag = 1;
n = size(S,1);

M_all = zeros(size(S));
[LS,~] = pca1(S,3);
for i = 1:length(idx)
%     [chuzhiA,~] = nndsvd(S,i,flag);
%     params.tol = 10^(-6);
%     params.Hinit = chuzhiA;
%     
%     [H_all,~,~] = symnmf_newton(S,i,params);
%     [~,indx] = max(H_all, [], 2);
    indx = kmeans(LS,i,'Start',zeros(i,3));
    M = consen(indx);
    M_all = M_all + M;
end

M_all(M_all <= tau*length(idx)) = 0;
M_all = (1/2)*(M_all + M_all');

D = diag(M_all*ones(n,1));
Prw = eye(n) - D^(-1/2)*M_all*D^(-1/2);
% Prw = eye(n) - D^(-1)*M_all;
% P = (eye(n)/D)*M_all;
% all_eigs = eig(P);
if n>=1000
    No_eigs = 100;
    all_eigs = real(eigs(Prw,No_eigs,'sm'));
else
    all_eigs = real(eig(Prw));
end

    function M = consen(indx)
        
        % [m,~] = size(S);
        m = length(indx);
        
        M = zeros(m);
        
        for ii = 1:m
            for j = 1:m
                if indx(ii)==indx(j)
                    M(ii,j) = 1;
                end
            end
        end
    end
end

function [T] = sc_dpg(X, g, y)

%X = log(1+sc_norm(X));

[setmatrx, setnames, setgenes] = pkg.e_getgenesets;
[~,ix,iy]=intersect(setgenes,g);

setmatrx=setmatrx(:,ix);   % s x g
%setgenes=setgenes(ix);
%g=g(iy);

%isequal(gcm,g)
%isequal(gcm,setgenes)
X=X(iy,:);        % g x c

Z=setmatrx*X;   % s x c

p_val=ones(size(Z,1),1);
avg_log2FC=nan(size(Z,1),1);
for k=1:size(Z,1)
    a=Z(k,y);    
    b=Z(k,~y);
    p_val(k) = ranksum(a,b);
    if ~isnan(p_val(k))
        [ax]=nbinfit(a);
        [bx]=nbinfit(b);
        %[ax(1) bx(1)]
        avg_log2FC(k) = ax(1)./bx(1);    
    end
end
if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
end
T=table(setnames, avg_log2FC, p_val, p_val_adj);
T(isnan(T.p_val),:)=[];
T=sortrows(T,4);

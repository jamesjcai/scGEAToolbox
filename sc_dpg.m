function [T] = sc_dpg(X, Y, g, setmatrx, setnames, setgenes)

%X = log(1+sc_norm(X));
if nargin<6
    [setmatrx, setnames, setgenes] = pkg.e_getgenesets;
end
[~,ix,iy] = intersect(upper(setgenes),upper(g));

setgenes=setgenes(ix);
setmatrx=setmatrx(:,ix);    % s x g

X=X(iy,:);              % g x c
Y=Y(iy,:);              % g x c
g=g(iy);

Zx = setmatrx*X;               % s x c
Zy = setmatrx*Y;               % s x c

gsetsize = sum(setmatrx,2);   % gene number per set

n=size(Zx,1);

p_val = ones(n,1);
avg_log2FC = nan(n,1);
v1 = nan(n,1);
v2 = nan(n,1);
n1 = nan(n,1);
n2 = nan(n,1);
m1 = nan(n,1);
m2 = nan(n,1);


for k=1:n
    if any(setmatrx(k,:))
    a=Zx(k,:);
    b=Zy(k,:);
    p_val(k) = ranksum(a,b);
    if ~isnan(p_val(k)) && p_val(k)<1e-3
        %[ax]=nbinfit(a);
        %[bx]=nbinfit(b);
        [ax]=mean(a);
        [bx]=mean(b);
        avg_log2FC(k) = log2(ax(1)./bx(1));
        v1(k)=ax(1);
        v2(k)=bx(1);
        n1(k)=numel(a);
        n2(k)=numel(b);
        m1(k)=sum(a>0);
        m2(k)=sum(b>0);
    end
    end
end
%warning on
if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
end

T=table(setnames, gsetsize, v1, v2, avg_log2FC, m1, n1, m2, n2, p_val, p_val_adj);
T(isnan(T.p_val)|isnan(T.avg_log2FC)|abs(T.avg_log2FC)<1,:)=[];
T = sortrows(T, 'p_val_adj', 'ascend');
T=T(T.p_val_adj<0.01 & T.gsetsize>=5,:);

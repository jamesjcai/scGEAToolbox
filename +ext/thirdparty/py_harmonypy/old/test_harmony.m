t1=readtable('pcs.csv');
s=table2array(t1);
N=size(s,1);

t2=readtable('meta.csv');
c_batch_id=string(t2.dataset);
[c,cL]=grp2idx(c_batch_id);
phi=dummyvar(c)';
phi_n=numel(unique(c));
N_b = sum(phi,2);

Pr_b = N_b / N;

%%
figure;
scatter3(s(:,1),s(:,2),s(:,3),[],c)
colormap(lines(phi_n))

%%
t3=readtable('adj.tsv','FileType','text');
s2=table2array(t3)';
figure;
scatter3(s2(:,1),s2(:,2),s2(:,3),[],c)
colormap(lines(phi_n))

%%

data=readtable("pcs.csv");
s=table2array(data);

Z_cos=s./max(s,[],2);
Z_cos=Z_cos./vecnorm(Z_cos,2,2);

Z_corr=s;
Z_orig=s;


N=size(s,1);
K=min([round(N/30) 100]);
B=size(phi,1);  % # number of batch variables
d=size(Z_cos,2);

[idx,C,sumd]=kmeans(Z_cos,K,'MaxIter',10);

Y=C';
Y=Y./vecnorm(Y,2,2);
dist_mat=2*(1-Y'*Z_cos');
sigmav=0.1*ones(K,1);

R = -dist_mat;
R = R ./ sigmav;
R=R-max(R);
R=exp(R);
R = R ./sum(R);
E=sum(R,2)*Pr_b';
O=R*phi';
E=sum(R,2)*Pr_b';

% O=R'*phi;
% https://github.com/slowkow/harmonypy/blob/master/harmonypy/harmony.py


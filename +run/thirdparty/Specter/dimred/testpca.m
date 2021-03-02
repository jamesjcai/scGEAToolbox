n = 100000;
m = 200;
density = 0.2;
k = 5;
Y = sprand(n, m, density);
Y = 100*full(Y);

%  tic;
%  [c,s,lambda] = pca(Y);
%  toc;
%  
tic;
[u,v] = fsvd(Y,k,3);
toc;
u(1:5,1:5), v(1:5, 1:5)
tic;
[u, v] = svds(Y,k);
toc;
/* u(1:5,1:5), v(1:5, 1:5) */
/* tic; */
/* [u,v,s] = svdsecon(Y,k); */
/* toc; */
/* u(1:5,1:5), v(1:5, 1:5) */

opts.tol = 1e-8;
opts.maxit = 150;


%  s(1:10,1:10), score(1:10,1:10)

%tic;
%maxIter = 1000;
%numRep = 20;
%litekmeans(Y,10,'MaxIter', maxIter, 'Replicates', numRep);
%toc;

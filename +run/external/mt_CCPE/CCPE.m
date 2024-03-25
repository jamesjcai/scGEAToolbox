function [pseudotime] = CCPE(X,lambda,gamma,sigma)
% [mappedX, mapping] = compute_mapping(X', 'PCA', 3);
% W=mapping.M;
% Z=mappedX;
[W,Z] = pca(X','NumComponents',3);
Z=Z';
N=size(Z,2);
k=N;
vaxis='x';
%initial Y
Y=Z;

iter = 0;
while iter < 200
      [~, ~,Z_p,~,~,~] = HelixFit(Z,vaxis);
      [R] = get_R(Z,Y,sigma);
      T=diag(diag(ones(k,N)*R));
      Q=inv((1+lambda+gamma)*eye(N,N)-gamma*R*inv(T)*R');
      A=Z_p*Q*X';
      [U,~,V] = svd(A);
      I=eye(size(W));
      W=V*I*U';
      Z=(W'*X+lambda*Z_p)*Q;
      Y=Z*R*inv(T);
      error1=trace((X-W*Z)*(X-W*Z)');
      error2=lambda*trace((Z-Z_p)*(Z-Z_p)');
      error3=gamma*(trace(Z*Z')-2*trace(R'*Z'*Y)+trace(Y*T*Y'));
      MSE_error=error1+error2+error3;
      iter=iter+1;
      if iter > 199
         fprintf('iter achieved!\n');
         break;
      end     
end

pseudotime=Z(1,:)';
end


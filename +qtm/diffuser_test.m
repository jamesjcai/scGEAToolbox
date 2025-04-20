N = 4;
diffuserG = qtm.diffuser(N);
plot(diffuserG)



e = ones(2^N,1);
e = e/norm(e);
M = eye(2^N) - 2*(e*e');   % reflection matrix
norm(getMatrix(diffuserG) - M)


 function [Y,w, L_1,L_2,L_3] = Algorithm_4(x_1, x_2, x_3, sample_num, iter_num, lam_1, d_num)
% x_1, x_2 and x_3 is the local sample-spectrum
% Y is the global sample-spectrum

w=ones(3*sample_num,1)./(3*sample_num);

Y_final = [];
w_final = [];

final_err= inf;

%% doing Algorithm 3
for iter=1:iter_num

w_1 = diag(sqrt(w(1:sample_num,1)));
w_2 = diag(sqrt(w((sample_num+1):2*sample_num,1)));
w_3 = diag(sqrt(w((2*sample_num+1):3*sample_num,1)));


s_1 = sum(diag(((eye(sample_num)-(x_1)\(x_1)))*((eye(sample_num)-(x_1)\(x_1)))'));
s_2 = sum(diag(((eye(sample_num)-(x_2)\(x_2)))*((eye(sample_num)-(x_2)\(x_2)))'));
s_3 = sum(diag(((eye(sample_num)-(x_3)\(x_3)))*((eye(sample_num)-(x_3)\(x_3)))'));

M = (1/s_1)*(w_1*(eye(sample_num)-(x_1*w_1)\(x_1*w_1)))*(w_1*(eye(sample_num)-(x_1*w_1)\(x_1*w_1)))'+(1/s_2)*(w_2*(eye(sample_num)-(x_2*w_2)\(x_2*w_2)))*(w_2*(eye(sample_num)-(x_2*w_2)\(x_2*w_2)))'+(1/s_3)*(w_3*(eye(sample_num)-(x_3*w_3)\(x_3*w_3)))*(w_3*(eye(sample_num)-(x_3*w_3)\(x_3*w_3)))';

% compute the global sample-spectrum (Y) based on the eigenvalue decomposition
[Y,Eigen_Value_all]=Find_K_Min_Eigen(M,d_num+1);
Y=Y(:,2:(d_num+1))';

L_1 = (Y*w_1)/(x_1*w_1);
L_2 = (Y*w_2)/(x_2*w_2);
L_3 = (Y*w_3)/(x_3*w_3);

err = sum( diag(Y*M*Y') );
if err-final_err<0
    final_err = err;
    Y_final = Y;
    w_final = w;
else
    break
        
end
%% obtain err
[ E_1,E_2,E_3 ] = compute_err( w, sample_num, Y, L_1, L_2, L_3, x_1, x_2, x_3, s_1, s_2, s_3 );

%% doing Algorithm 2, solve w
w_1= Algorithm_2( sample_num, w, E_1, E_2, E_3, lam_1);
w=w_1;
 


end

Y=Y_final;
w=w_final;


 end



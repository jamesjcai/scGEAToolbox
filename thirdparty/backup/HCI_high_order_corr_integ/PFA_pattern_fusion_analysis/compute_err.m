function [ E_1,E_2,E_3 ] = compute_err( w, sample_num, Y, L_1, L_2, L_3, x_1, x_2, x_3, s_1, s_2, s_3  )

E_1 = [];
E_2 = [];
E_3 = [];

for i=1:sample_num
    E_1 =[E_1; (1/s_1)*w(i,1)*(Y(:,i)-L_1*x_1(:,i))'*(Y(:,i)-L_1*x_1(:,i))];
    E_2 =[E_2; (1/s_2)*w(sample_num+i,1)*(Y(:,i)-L_2*x_2(:,i))'*(Y(:,i)-L_2*x_2(:,i))];
    E_3 =[E_3; (1/s_3)*w(2*sample_num+i,1)*(Y(:,i)-L_3*x_3(:,i))'*(Y(:,i)-L_3*x_3(:,i))];
 
end

end


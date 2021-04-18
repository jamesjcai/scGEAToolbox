pw1 = fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'PFA_pattern_fusion_analysis');
addpath(pth);

load ../../example_data/example10xdata.mat
%%
%[X]=run_magic(X,true);
%[X]=sc_norm(X,'type','deseq');
%celllist=string(1:size(X,2));
%%
F1=corr(X);
F2=corr(F1);



%% capture the local sample-spectrum for each biological data type by Algorithm_1
data_1 = X;
data_2 = F1;
data_3 = F2;

sample_num = size(data_1,2);

[u_1,eig_v_1]=Algorithm_1(data_1,sample_num);
[u_2,eig_v_2]=Algorithm_1(data_2,sample_num);
[u_3,eig_v_3]=Algorithm_1(data_3,sample_num);

d_1 =1;
for d_1=1:sample_num
    if sum(eig_v_1(1:d_1))/sum(eig_v_1)>0.8
        break;
    end
end

d_2 =1;
for d_2=1:sample_num
    if sum(eig_v_2(1:d_2))/sum(eig_v_2)>0.9
        break;
    end
end

d_3 =1;
for d_3=1:sample_num
    if sum(eig_v_3(1:d_3))/sum(eig_v_3)>0.9
        break;
    end
end
  
[u_1,eig_v_1]=Algorithm_1(data_1,d_1);
[u_2,eig_v_2]=Algorithm_1(data_2,d_2);
[u_3,eig_v_3]=Algorithm_1(data_3,d_3);

x_1 = u_1'*data_1; %% the local sample-spectrum
x_2 = u_2'*data_2; %% the local sample-spectrum
x_3 = u_3'*data_3; %% the local sample-spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% capture the global sample-spectrum according to Algorithm_4
d_num = min([d_1,d_2,d_3]);

iter_num =1000;
lam_1 = 1;

[Y,w,L_1,L_2,L_3] = Algorithm_4(x_1, x_2, x_3, sample_num, iter_num,  lam_1, d_num);%% Y is the global sample-spectrum

% xlswrite('res_fileNameaa',Y);


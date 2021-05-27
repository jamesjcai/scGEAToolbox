% Estimation of beta for the k_NN, for several k and dimensions and gamma=1

gamma = 1;

k_min = 1;
k_max = 7;

d_min = 2;
d_max = 7;

N = 10000;
n_trials = 50;


beta_hat = zeros(d_max-d_min+1,k_max-k_min+1);  % array of estimates of beta
beta_std = beta_hat;                            % (estimated) standard deviation of beta estimates

tic

for d = d_min:d_max
   %clear X Xr
   display(['Dimension: ',num2str(d)]);
   
   for k = k_min:k_max
       display(['k: ',num2str(k)])
       knnlen1_norm = 0; knnlen2_norm = 0;
       
       for i = 1:n_trials
           X = rand(d,N);
           Xr = X(:);
           [L Graph] = kNNgraphmex(Xr,N,d,k,gamma);
           
           %L/(N^(1-gamma/d))
           
           knnlen1_norm=knnlen1_norm+L/(N^(1-gamma/d));
           knnlen2_norm=knnlen2_norm+(L/(N^(1-gamma/d)))^2;
           
       end
       beta_hat(d-d_min+1,k-k_min+1) = knnlen1_norm/n_trials;
       beta_std(d-d_min+1,k-k_min+1) = sqrt((knnlen2_norm-(knnlen1_norm/n_trials)^2*n_trials)/(n_trials-1));
       
   end
   
end

toc


%errorbar(beta_hat,beta_std)
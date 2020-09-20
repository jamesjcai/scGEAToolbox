function [d,H,knnlenavg_vec,knnlenstd_vec] = knn_graph_estim_2(Dist, kneighbors, gamma, M, N, samp_points)

% Dimension and alpha-entropy estimatation of a data set using the k-NN graph
%
% Inputs: 
%   Dist = matrix of dissimilarities between points
%   kneighbors = number of neighbors in k-NN graph
%   gamma = edge weighting factor for k-NN
%   M = number of independent LS runs
%   N = number of resampling trials per LS dwell
%       M and N provide a bias-variance tradeoff for the estimates;
%       values of interest are of the form (M,N)=(M,1) or (M,N) = (1,N)
%   samp_points = vector of indices of LS dwells to be used 
% Outputs:
%   d = LS estimate of intrinsic dimension
%   H = LS estimate of entropy
%   knnlenavg_vec = avg length of k-NN over samp range
%   knnlenstd_vec = std of length of k-NN over samp range
% Required subroutines:
%   kNNlengthmex.c
%   beta_k_NN_estimate.mat (attention, only has beta constants for gamma = 1!)
%
%
% Written by Jose Costa (jcosta@umich.edu), Alfred Hero (hero@eecs.umich.edu), 2004
%
%
% BEGIN Copyright notice
% 
% Matab scripts for intrinsic dimension and entropy estimation using k-nearest
% neighbor graphs. The details of the algorithms can be found in:
% 
%   J. A. Costa and A. O Hero, "Entropic Graphs for Manifold Learning",
%   Proc. of IEEE Asilomar Conf. on Signals, Systems, and Computers,
%   Pacific Groove, CA, November, 2003.
%
%   J. A. Costa and A. O. Hero, "Geodesic Entropic Graphs for Dimension and
%   Entropy Estimation in Manifold Learning", 
%   to appear in IEEE Trans. on Signal Processing, Aug., 2004. 
%
% Published reports of research using the code provided here (or a modified version)
% should cite the two articles referenced above.
%
% Comments and questions are welcome. We would also appreciate hearing about 
% how you used this code, improvements made to it, etc. You are free to modify the
% code, as long as you reference the original contributors.
%
% END Copyright notice


clear beta_hat beta_std
load beta_k_NN_estimate

tsamp=size(Dist,1);

Q = length(samp_points);
knnlenavg_vec = zeros(M,Q);
knnlenstd_vec = zeros(M,Q);

for i = 1:M
    
    % Perform resampling estimation of mean knn length
    
    j = 1;
    for n = samp_points
        knnlen1=0;  knnlen2=0;
	    for trial = 1:N;
		    indices = randperm(tsamp);
		    indices = indices(1:n);
            Dr=Dist(indices,:);
            Drr=Dr(:,indices);
            L = knn_distmat(Drr,kneighbors);
            knnlen1 = knnlen1+sum(L);
            knnlen2 = knnlen2+sum(L)^2;
	    end
        knnlenavg_vec(i,j) = knnlen1/N;
        if (N ~= 1) %
            knnlenstd_vec(i,j)=sqrt((knnlen2-(knnlen1/N)^2*N)/(N-1));
        else %
            knnlenstd_vec(i,j)=0; %
        end %
        j = j+1;
    end


    %Compute LS estimate of d and H
    
    %Q = length(samp_points);
    A = [log(samp_points)',ones(Q,1)];

    sol1 = inv(A'*A)*A'*log(knnlenavg_vec(i,:))';
    %d=round(gamma/(1-sol1(1)));
    dvec(i) = gamma/(1-sol1(1));

    if ((kneighbors <= size(beta_hat,2)) & (round(dvec(i)) <= size(beta_hat,1)+1) & (round(dvec(i)) > 1) )
        % assumes gamma = 1
        %H=((sol1(2)-log(beta_hat(d-1,kneighbors)))/(1-sol1(1)))*log2(exp(1)); % express in terms of bits
        Hvec(i) = (dvec(i)/gamma)*(sol1(2)-log(beta_hat(round(dvec(i))-1,kneighbors)))*log2(exp(1)); % express in terms of bits
    else
        display(['Beta(m = ',num2str(dvec(i)),', k = ',num2str(kneighbors),') not available for entropy estimation'])
        Hvec(i) = -10;
    end
    
    d = mean(dvec);
    H = mean(Hvec);  
    
end
function M = selective_sampling(clusters, k, L)
 % fprintf('------------------------------- selective sampling -------------------------------------\n');  
   [m, n] = size(clusters);
   M = []; % store indices
   for iter =1:1:L
       for i = 1:1:m 
           for j=1:1:k 
               I = find(clusters(i,:)==j);
               s = I(randperm(length(I), 1));
               if (~ismember(s,M))
                   M = [M s];
               end
               if length(M) >= L 
                   return;
               end
           end
       end 
   end
   if (length(M) < L)
       C = setdiff(1:n, M);
       sC = C(1:(L-length(M)));
       M = [M sC];
   end
end

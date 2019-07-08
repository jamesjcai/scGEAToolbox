% function cl = clkmeans(x,k,sfct,oldcl)
%
% DESCRIPTION
%   provides cluster labels 1 to k applying k-means to x
%   make sure samples (rows) in x are randomly permuted since the
%   first k samples (rows) are used for initialization
%   if oldcl is given their centroids are used to initalize centroids
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = clkmeans(x,k,sfct,oldcl)

if exist('oldcl')
   if (length(oldcl)==size(x,1))&(max(oldcl)==k),
      disp('clkmeans: initializing means with centroids of given clustering')
      [centers options cln] = kmeans(clucent(x,oldcl),x,zeros(1,14),sfct);
   else
      r = randperm(size(x,1));
      initcenters = r(1:k);
      disp(['clkmeans: initializing means with random samples ' num2str(initcenters)]);
      [centers options cln] = kmeans(x(initcenters,:),x,zeros(1,14),sfct);
   end;
else
   disp('clkmeans: initializing means with k first samples'); 
   [centers options cln] = kmeans(x(1:k,:),x,zeros(1,14),sfct);
end;
[x y] = find(cln);
cl = sortrows([x y]);
cl = checkcl(cl(:,2)');


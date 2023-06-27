% function [centres, options, post, errlog] = kmeans(...)
%
% This function has been adopted from a k means implementation
% copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

function [centres, options, post, errlog] = XX_kmeans(centres, data, options, measure)

[ndata, data_dim] = size(data);
[ncentres, dim] = size(centres);

if dim ~= data_dim
  error('Data dimension does not match dimension of centres')
end

if (ncentres > ndata)
  error('More centres than data')
end

% Sort out the options
if (options(14))
  niters = options(14);
else
  % niters = 25; % AS was this for a long time
  niters = 50;
  niters = 100; % just for ml paper examples
end

store = 0;
if (nargout > 3)
  store = 1;
  errlog = zeros(1, niters);
end

% Check if centres and posteriors need to be initialised from data
if (options(5) == 1)
  % Do the initialisation
  perm = randperm(ndata);
  perm = perm(1:ncentres);

  % Assign first ncentres (permuted) data points as centres
  centres = data(perm, :);
end
% Matrix to make unit vectors easy to construct
id = eye(ncentres);

% Main loop of algorithm
for n = 1:niters

  % Save old centres to check for termination
  old_centres = centres;
  
  % Calculate posteriors based on existing centres
  d2 = sqrt(-log((feval(measure,data,centres))));
  % Assign each point to nearest centre
  [minvals, index] = min(d2', [], 1);
  post = id(index,:);

  num_points = sum(post, 1);
  % Adjust the centres based on new posteriors
  for j = 1:ncentres
    if (num_points(j) > 0)
      centres(j,:) = sum(data(find(post(:,j)),:), 1)/num_points(j);
    end
  end

  % Error value is total squared distance from cluster centres
  e = sum(minvals);
  if store
    errlog(n) = e;
  end
  if options(1) > 0
    fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
  end

  if n > 1
    % Test for termination
    if max(max(abs(centres - old_centres))) < options(2) & ...
        abs(old_e - e) < options(3)
      options(8) = e;
      return;
    end
  end
  old_e = e;
end

% If we get here, then we haven't terminated in the given number of 
% iterations.
options(8) = e;
if (options(1) >= 0)
%  disp('Warning: Maximum number of iterations has been exceeded');
end


% function x = norml(x,mode)
%
% copyright (c) 1998-2011 by Alexander Strehl

function x = norml(x,mode)

switch mode
  case 0,
  case 1,
   x = diag(sparse(1./sum(abs(x),2))) * x;
  case 2,
   x = diag(sparse(1./sum((x.^2),2).^(1/2))) * x;
  otherwise, 
     disp('norml: mode not supported');
end; 
 

% function cl = clrand(x,k,sfct)
%
% DESCRIPTION
%   provides cluster labels 1 to k drawn from a uniform distribution
%   however it makes sure each label occurs at least once
%   sfct is ignored, x is only used to figure the number of objects 
%
% Copyright (c) 1998-2011 by Alexander Strehl

function cl = clrand(x,k,sfct)

if k<=size(x,1)
  cl = floor(k*rand(1,size(x,1)))+1;
  if (max(onetomax(cl)) ~= k)
    rindex = randperm(size(x,1));
    cl(rindex(1:k)) = 1:k;
  end;
else
  disp('clrand: more clusters than samples requested');
  cl = 1:size(x,1);
end;

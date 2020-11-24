% function imageshq(s)
%
% copyright (c) 1998-2011 by Alexander Strehl

function imageshq(s)

imagesc(s);
map = colormap;

if (sum(sum((map(:,1)~=map(:,2))+(map(:,1)~=map(:,3)))==0)),
   
   n = 100;
   if (size(s,1)>n),
      subsample = randperm(size(s,1));
      subsample = subsample(1:n);
      s = s(subsample,subsample);
      disp(['imageshq-warning: using subsample of size ' num2str(n) ' to determine transformation']);
   end;
   
   map = gray;
   cols = size(map,1);
   range = 0:1/(cols-1):1;
   
   c = hist(reshape(s,prod(size(s)),1),range);
   ca = cumsum(c);
   badones = find(ca==0);
   ca(badones) = ones(size(badones));
   
   cb = prod(size(s))*range;
   badones = find(cb==0);
   cb(badones) = ones(size(badones));
   
   newmap = map./((ones(3,1)*(cb./ca))');
   colormap(newmap);
   
   disp('imageshq-warning: colormap has been modified for histeq');
   
else
   disp('imageshq-error: not a graylevel map');
end;

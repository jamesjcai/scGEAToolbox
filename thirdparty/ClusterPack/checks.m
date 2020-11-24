% function s = checks(s)
% 
% DESCRIPTION
%   Checks a similarity matrix for validity and fixes
%   detected problems
%
% Copyright (c) 1998-2011 by Alexander Strehl


function s = checks(s)
if ~isempty(s),
   if sum(sum(~isreal(s))),
      disp('checks: complex similarities found');
      s = real(s);
      disp('checks: using real component');
   end;
   if size(s,1)~=size(s,2)
      disp('checks: s is not square');
      s = s(1:min(size(s)),1:min(size(s)));
      disp('checks: using indrawn square');
   end;
   mas = max(max(s));
   mis = min(min(s));
   if (mas>1)|(mis<0),
      disp(['checks: similarity more than 1 or less than 0 detected: values ' num2str(mis) ' to ' num2str(mas) ]); 
      s(find(s<0))=0;
      s(find(s>1))=1;
      disp('checks: bounded');
   end;
   if sum(sum(isinf(s)|isnan(s))),
      disp('checks: non-finite similarity detected !!! (serious)'); 
      if 0,
         s(find(isinf(s)|isnan(s))) = 0; % hangs the computer - no idea why...
      else
        [a,b] = find(isfinite(s));
         c = find(isfinite(s));
         s = sparse(a,b,s(c)); 
      end;
      disp('checks: made zero !!! (serious)');
   end;
   if (sum(sum(s~=s'))>0),
      disp('checks: s is not symmetric');
      s = (s+s')./2;     
      disp('checks: symmetrised');
   end;
   if (size(s,1)==size(s,2)),
      if sum(diag(s)~=1)
         disp('checks: self-similarity s(i,i) not always 1');
         for i=1:size(s,1),
            s(i,i) = 1;
         end;
         disp('checks: diagonal made 1');
      end;
   end;
else
   disp('checks: empty similarity matrix');	
end;

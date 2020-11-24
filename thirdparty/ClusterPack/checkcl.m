% function cl = checkcl(cl)
%
% DESCRIPTION
%   Checks a cluster labeling for validity and fixes
%   detected problems
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = checkcl(cl)

if ~isempty(cl),
   if min(cl)<=0,
      disp('checkcl: non-positive values in cl');
      cl = cl-min(cl)+1;
      disp('checkcl: offset to minimum 1');   
   end;
   if max(cl) ~= max(onetomax(cl)),
      disp('checkcl: not a dense integer mapping');
      cl = onetomax(cl);
      disp(['checkcl: remapped to 1 to ' num2str(max(cl))]);
   end;
else
   disp('checkcl: empty clustering');
end;

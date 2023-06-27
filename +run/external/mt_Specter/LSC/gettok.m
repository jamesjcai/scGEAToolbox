% Author: Modified package by Van Hoan Do
function [word, remainder] = gettok(string)
%
%  function [word, remainder] = gettok(string)
%
%  Retrieves the first blank separated token from the string.
%
si = findstr(string,' ');
lstring=length(string);
lsi=length(si);
if ( lsi == 0 )
   word=string;
   remainder='';
   return
end
firstb=si(1);
if ( firstb > 1 )
  word=string(1:firstb-1);
  remainder=string(firstb+1:lstring);
  return;
end
tmp=1;
while ( tmp < lsi )
  if ( si(tmp+1) == si(tmp)+1 )
     tmp=tmp+1;
  else
     break;
  end
end
if ( tmp == lstring )
  word=-1;
  remainder=string;
  return;
end
word=string(si(tmp)+1:si(tmp+1)-1);
remainder=string(si(tmp+1)+1:lstring);

    
   

% function p = evalmse(trueclass,cl,x,sfct)
% 
% DESCRIPTION
%   quality based on average error per object and per feature
%   ignores trueclass - just for compatibility
%
% Copyright (c) 1998-2011 by Alexander Strehl

function p = evalmse(trueclass,cl,x,sfct)

c = clucent(x,cl);
sse = 0;
for i = 1:max(cl),
   ins = find(cl==i);
   e = x(ins,:)-(ones(length(ins),1)*c(i,:));
   sse = sse + sum(sum(e.*e));
end;
p = exp(-sse/prod(size(x)));

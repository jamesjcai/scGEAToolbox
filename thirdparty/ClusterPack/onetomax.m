% function c2 = onetomax(c1)
%
% copyright (c) 1998-2011 by Alexander Strehl

function c2 = onetomax(c1)

temp = sortrows(mapdense([c1' (1:length(c1))'],1),2);
c2 = temp(:,1)';



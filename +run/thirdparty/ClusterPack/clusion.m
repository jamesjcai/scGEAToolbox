% function h = clusion(s,cl,hieq,sq)   
%
% DESCRIPTION
%   - displays the CLUSION cluster visualization plot of a similarity matrix
%   - returns a column vector of handles h to the cluster separating lines
% ARGUMENTS
%   s    - square similarity matrix (n x n)
%   cl   - cluster label row vector (1 x n)
%          entries must be integers from 1 to k 
%   hieq - 0/1 to indicate if histogram should be equalized 
%          (optional, default 0) 
%   sq   - 0/1 to indicate if figure axis should be square
%          (optional, default 0)
% EXAMPLE
%   clusion;
%   clusion([1 0 .5; 0 1 0; .5 0 1],[1 2 1]);
% REFERENCE
%   please refer to the following paper if you use CLUSION:
%     A. Strehl and J. Ghosh, "Relationship-based Clustering and
%     Visualization for High-dimensional Data Mining", version 2.0, 2011/05/01
%     Issue on Mining Web-based Data for E-Business Applications
%     of the INFORMS Journal on Computing, 2002
% RELEASE
%   version 2.0, 2011/05/01, tested on WIN7 Octave 3.2.4 and LNX86 Matlab 5.2.0.3084
%   available from http://www.strehl.com
%   license granted for research use ONLY (see README)
%   copyright (c) 1998-2011 by Alexander Strehl

function h = clusion(s,cl,hieq,sq)

if exist('histeq.m','file')
   profihisteq = 1;
else
   profihisteq = 0;
end

if ~exist('s')
   disp('clusion-warning: no arguments - displaying illustrative example:');
   s=[1 0.0101444233725719 0.666666666666667 0.000324894409316972 0.000144423583260196 0 0.000162460400277432;0.0101444233725719 1 0.0204877164355 0.00700161873707528 0.00504709317581638 0 0.00728303954110406;0.666666666666667 0.0204877164355 1 0.000649366867304378 0.00028876375792712 0 0.000324894409316972;0.000324894409316972 0.00700161873707528 0.000649366867304378 1 9.99900009999e-005 0 0.666666666666667;0.000144423583260196 0.00504709317581638 0.00028876375792712 9.99900009999e-005 1 0.99985557641674 0.000103983778530549;0 0 0 0 0.99985557641674 1 0;0.000162460400277432 0.00728303954110406 0.000324894409316972 0.666666666666667 0.000103983778530549 0 1]
   cl=[3 2 3 2 1 1 2]
   disp('clusion-advice: type "help clusion" for information about usage');
end

if size(s,1)~=size(s,2)
   disp('clusion-error: similarity matrix not square');
   h = [];
   return;
end

if ~exist('cl')
   disp('clusion-warning: no clustering given');
   cl = ones(1,size(s,1));
end

if size(cl,1)~=1
   disp('clusion-error: clustering must be row vector');
   return;
end

if size(s,1)~=length(cl)
   disp('clusion-error: clustering and similarity matrices mismatch in size');
   return;
end

if ~exist('hieq')
   hieq = 0;
end

if ~exist('sq')
   sq = 0;
end

colormap gray;

index = [1:length(cl); cl];
index = sortrows(index',2)';
order = index(1,:);
ticks = conv(index(2,:),[1,-1]); 
ticks(1) = 0; ticks(size(index,2)+1) = 0;
where = [1 find(ticks) length(cl)+1]-.5;
xline = [where;where];
yline = [ones(1,length(where))-.5; (length(s)+.5)*ones(1,length(where))];

if (hieq==0)
   imagesc(s(order,order));
else
   if profihisteq
      imagesc(histeq(full(s(order,order))));
   else
      imageshq(s(order,order));
   end
end
handa = line(xline,yline);
handb = line(yline,xline);
set([handa handb],'color','yellow');
set([handa(1) handa(length(handa)) handb(1) handb(length(handb))],'color','black');


if (length(cl)<=25)
   renumtick(gca,cl,cl);
else
   axis off;
end



if sq
   axis square;
else
   axis normal;
end

h = [handa;handb];

colormap(1-colormap);
set(h,'Color',[1.0 0.1 0.1]);


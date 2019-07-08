% function [cl, index] = clbtocl(clb,p)
%
% copyright (c) 1998-2011 by Alexander Strehl

function [cl, index] = clbtocl(clb,p)
randbreakties = 1;

allzerocolumns = find(sum(clb,1)==0);
if ~isempty(allzerocolumns),
   disp(['clbtocl: ' num2str(length(allzerocolumns)) ' objects (' num2str(100*length(allzerocolumns)/size(clb,2),'%.0f') '%) with all zero associations']);
   clb(:,allzerocolumns) = rand(size(clb,1),length(allzerocolumns));
end;
if randbreakties,
   clb = clb + rand(size(clb))/10000;
end;

clb = norml(clb',1)';
m = max(clb,[],1);
cl = zeros(1,size(clb,2));
winnersprop = zeros(1,size(clb,2));
for i=size(clb,1):-1:1,
  a = find(m==clb(i,:));
  cl(a) = i*ones(1,length(a));
  winnersprop(a) = clb(i,a);
end;

if ~exist('p')
  index = []; 
else
  index = find(winnersprop<p); 
end;
disp(['clbtocl: delivering ' num2str(max(cl)) ' clusters']);
disp(['clbtocl: average posterior prob is ' num2str(mean(winnersprop)) ]);
if (length(cl)<=7),
  disp('clbtocl: winning posterior probs are');
  disp(winnersprop);
  disp('clbtocl: full posterior probs are');
  disp(clb);
end;


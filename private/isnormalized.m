function [y]=isnormalized(X)
%This function assumes that unnormalized matrix contains all integers.
y=~all(mod(X,1)==0,'all');
% y = all(X==floor(X),'all')
% [e,ker] = ee_error(X,Wp,Wn,l) EE objective function error
function [e,ker] = ee_error(X,Wp,Wn,l)
sqd = sqdist(X); ker = exp(-sqd);
e = Wp(:)'*sqd(:) + l*Wn(:)'*ker(:);

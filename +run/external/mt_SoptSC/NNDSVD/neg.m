function [Am] = neg(A);
Am = (A<0).*(-A);
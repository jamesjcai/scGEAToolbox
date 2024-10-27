function [H, P] = Hbeta(D, beta)
D = (D - min(D)) / (max(D) - min(D) + eps);
P = exp(-D*beta);
sumP = sum(P);
H = log(sumP) + beta * sum(D.*P) / sumP;
P = P / sumP;
end
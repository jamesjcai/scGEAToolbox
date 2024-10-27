function [s] = ml_TSEE(X, t, dim)
if nargin < 3, dim = 2; end

% https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5477-8
% TSEE: an elastic embedding method to visualize the dynamic gene
% expression patterns of time series single-cell RNA sequencing data

assert(size(X, 2) == length(t))

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, 'external', 'ml_TSEE');
if ~(ismcc || isdeployed), addpath(pth); end

X = sc_transform(X);


Y = X';
y_time = t;


% Y - cell x gene matrix
% y_time - time

% input gene expression profile matrix Y, each row represents a cell, and each column represents a gene;
% Y = double(Y);
Y = Y - min(Y(:));
Y = Y ./ max(Y(:));

% input the time of cells as y_time;
y_time = double(y_time);
y_time = (y_time - min(y_time)) / (max(y_time) - min(y_time));

% calculate two weight matrices
N = size(Y, 1);
beta = 10;
Wp = ea(Y, 20);
Wp = (Wp + Wp') / 2;
Wp(1:N+1:N^2) = 0;
Wp = Wp / sum(Wp(:));

Wn = sqdist(Y) + beta * absdist(y_time);
Wn = (Wn + Wn') / 2;
Wn(1:N+1:N^2) = 0;
Wn = Wn / sum(Wn(:));

% initialize and perform TSEE
d = dim;
opts.maxit = 100;
opts.tol = 1e-3;
s = RandStream('mcg16807', 'Seed', 29);
RandStream.setGlobalStream(s);
l = 10;
[X, ~, ~, ~] = ee(Wp, Wn, d, l, opts); % X is the 2-dimensional embedding.

s = X;

end
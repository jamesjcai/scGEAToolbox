function [A] = pcrnet_denoised(X, varargin)
% Construct denoised GRN using tensor decomposition on bootstrapped PCR nets
%
% A = net.pcrnet_denoised(X);
% A = net.pcrnet_denoised(X, 'smplmethod', 'Jackknife');
%
% X is gene x cell matrix

import ten.*

if nargin < 1
    error(sprintf('USAGE: A=net.pcrnet_denoised(X);\n       A=net.pcrnet_denoised(X,''smplmethod'',''Jackknife'');'));
end

p = inputParser;
addOptional(p, 'smplmethod', "Jackknife", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["jackknife", "bootstrap"]));
addOptional(p, 'tdmethod', "CP", @(x) (isstring(x) | ischar(x)) & ismember(upper(string(x)), ["CP", "TUCKER"]));
addOptional(p, 'nsubsmpl', 10, @(x) fix(x) == x & x > 0);
addOptional(p, 'csubsmpl', 500, @(x) fix(x) == x & x > 0);
addOptional(p, 'savegrn', false, @islogical);
addOptional(p, 'donorm', true, @islogical);
parse(p, varargin{:});
tdmethod = p.Results.tdmethod;
nsubsmpl = p.Results.nsubsmpl;
csubsmpl = p.Results.csubsmpl;
smplmethod = p.Results.smplmethod;
savegrn = p.Results.savegrn;
donorm = p.Results.donorm;

switch upper(tdmethod)
    case "CP"
        tdmethod = 1;
    case "TUCKER"
        tdmethod = 2;
end
switch lower(smplmethod)
    case "jackknife"
        usebootstrp = false;
    case "bootstrap"
        usebootstrp = true;
end

if exist('@tensor/tensor.m', 'file') ~= 2
    error('Need Tensor Toolbox for MATLAB (https://www.tensortoolbox.org/)');
end

if exist('net.pcrnet', 'file') ~= 2
    error('Need net.pcrnet (scGEAToolbox)');
end

if donorm
    X = sc_norm(X, "type", "libsize");
    X = log1p(X);
end

[XM] = ten.i_nc(X, nsubsmpl, 3, csubsmpl, usebootstrp);

disp('Tensor decomposition')
[A] = ten.i_td1(XM, tdmethod);

if savegrn
    tstr = matlab.lang.makeValidName(string(datetime));
    save(sprintf('A_%s', tstr), 'A', '-v7.3');
end
end

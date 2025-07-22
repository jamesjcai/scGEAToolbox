function [Y] = e_locfit(X, t)

arguments
    X {mustBeNumeric,mustBeReal}
    t (:,1) {mustBeNumeric,mustBeReal}
end

%[t, idx] = sort(t);
%X = X(idx, :);

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'external', 'locfit');

% https://www.mathworks.com/matlabcentral/answers/3167-two-functions-with-the-same-name-how-to-directly-call-one-of-both
oldpth = pwd();
cd(pth);
predict = @predict;
cd(oldpth);

if ~(ismcc || isdeployed), addpath(pth); end
%pth = fullfile(pw1, '..', 'external', 'locfit', 'source');
%if ~(ismcc || isdeployed), addpath(pth); end
Y = zeros(size(X));
for k = 1:size(X, 2)
    fitm1 = locfit(t, X(:, k));
    Y(:, k) = predict(fitm1, t);
end
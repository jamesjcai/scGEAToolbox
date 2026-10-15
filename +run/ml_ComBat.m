function [X, Y] = ml_ComBat(X, Y)
pw1 = fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    combatcleanup = pkg.i_addpathtemp( ...
        fullfile(pw1, '..', 'external', 'ml_ComBat'));   %#ok<NASGU>
end

n1 = size(X, 2);
n2 = size(Y, 2);
batch = [ones(1, n1), 2 * ones(1, n2)];
mod = ones(size(batch))';
normxy = combat([X, Y], batch, mod);

X = normxy(:, 1:n1);
Y = normxy(:, n1+1:end);
end

function [scnout] = subnetwork_neighbors(obj, targetg, siz)
if nargin < 3, siz = 30; end
[~, idx] = ismember(targetg, obj.g);
w = 0.5 * (obj.A(idx, :) + obj.A(:, idx).');
[~, idxv] = maxk(abs(w), siz);
[scnout] = subnetwork(obj, [targetg; obj.g(idxv)]);
end
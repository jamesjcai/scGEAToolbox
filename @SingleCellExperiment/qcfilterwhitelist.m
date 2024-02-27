function obj = qcfilterwhitelist(obj, libszcutoff, mtratio, ...
    min_cells_nonzero, gnnumcutoff, whitelist)
if nargin < 6, whitelist = []; end
if nargin < 5 || isempty(gnnumcutoff), gnnumcutoff = 200; end
if nargin < 4 || isempty(min_cells_nonzero), min_cells_nonzero = 15; end
if nargin < 3 || isempty(mtratio), mtratio = 0.15; end
if nargin < 2 || isempty(libszcutoff), libszcutoff = 1000; end

if ~isempty(whitelist)
    assert(all(ismember(whitelist, obj.g)));
    [~, idxx] = ismember(whitelist, obj.g);
    Xresv = obj.X(idxx, :);
end

[~, keptg, keptidxv] = sc_qcfilter(obj.X, obj.g, libszcutoff, mtratio, ...
    min_cells_nonzero, gnnumcutoff);

for k = 1:length(keptidxv)
    obj = selectcells(obj, keptidxv{k});
end

[y] = ismember(obj.g, keptg);
obj.X = obj.X(y, :);
obj.g = obj.g(y);
[obj.X, obj.g] = sc_rmdugenes(obj.X, obj.g);

if ~isempty(whitelist)
    for k = 1:length(keptidxv)
        Xresv = Xresv(:, keptidxv{k});
    end
    [y, idxx] = ismember(whitelist, obj.g);
    
    if any(y)
        idxx=idxx(y);
        %assignin("base","X1",obj.X(idxx, :));
        %assignin("base","g1",obj.g(idxx));
        obj.X(idxx, :) = [];
        obj.g(idxx) = [];
    end
    %assignin("base","Xresv",Xresv)
    %assignin("base","whitelist",whitelist)
    obj.X = [obj.X; Xresv];
    obj.g = [obj.g; whitelist(:)];
end

end

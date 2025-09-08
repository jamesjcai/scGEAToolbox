function sce = qcfilterwhitelist(sce, libszcutoff, mtratio, ...
    min_cells_nonzero, gnnumcutoff, whitelist)
if nargin < 6, whitelist = []; end
if nargin < 5 || isempty(gnnumcutoff), gnnumcutoff = 200; end
if nargin < 4 || isempty(min_cells_nonzero), min_cells_nonzero = 15; end
if nargin < 3 || isempty(mtratio), mtratio = 0.15; end
if nargin < 2 || isempty(libszcutoff), libszcutoff = 1000; end

% if issparse(obj.X)
%     try
%         obj.X = full(obj.X);
%     catch
%     end
% end

if ~isempty(whitelist)
    assert(all(ismember(whitelist, sce.g)));
    [~, idxx] = ismember(whitelist, sce.g);
    Xresv = sce.X(idxx, :);
end

[~, keptg, keptidxv] = sc_qcfilter(sce.X, sce.g, libszcutoff, mtratio, ...
    min_cells_nonzero, gnnumcutoff);

for k = 1:length(keptidxv)
    sce = selectcells(sce, keptidxv{k}); % OK
end

[y] = ismember(sce.g, keptg);
sce.X = sce.X(y, :);
sce.g = sce.g(y);
[sce.X, sce.g] = sc_rmdugenes(sce.X, sce.g);

if ~isempty(whitelist)
    for k = 1:length(keptidxv)
        Xresv = Xresv(:, keptidxv{k});
    end
    [y, idxx] = ismember(whitelist, sce.g);
    
    if any(y)
        idxx=idxx(y);
        %assignin("base","X1",obj.X(idxx, :));
        %assignin("base","g1",obj.g(idxx));
        sce.X(idxx, :) = [];
        sce.g(idxx) = [];
    end
    %assignin("base","Xresv",Xresv)
    %assignin("base","whitelist",whitelist)
    sce.X = [sce.X; Xresv];
    sce.g = [sce.g; whitelist(:)];
end

end

function [X, genelist, keptidxv] = sc_qcfilter(X, genelist, libszcutoff, ...
    mtratio, min_cells_nonzero, gnnumcutoff)
%Basic QC filter
if nargin < 6, gnnumcutoff = 200; end
if nargin < 5 || isempty(min_cells_nonzero)
    min_cells_nonzero = 10; % 0.01
end
% if nargin<5 || isempty(dropout), dropout=1; end           % remove genes no expression
if nargin < 4 || isempty(mtratio), mtratio = 0.15; end % 0.10
if nargin < 3 || isempty(libszcutoff), libszcutoff = 500; end % 1000

try
    [X, keptidx] = sc_rmmtcells(X, genelist, mtratio);

catch ME
    warning(ME.message);
    Xobj = pkg.refwrap(X);
    clear X
    [keptidx] = obj_rmmtcells(Xobj, genelist, mtratio);
    X = Xobj.data;

end

keptidxv{1} = keptidx;
%if removemtgenes
%    [X,genelist]=sc_rmmtgenes(X,genelist);
%end
oldsz = 0;
newsz = 1;
c = 1;
while ~isequal(oldsz, newsz)
    oldsz = size(X);
    [X, genelist] = sc_filterg(X, genelist); % remove empty genes
    [X, keptidx] = sc_filterc(X); % remove empty cells
    keptidxv{end+1} = keptidx;
    [X, genelist] = sc_selectg(X, genelist, min_cells_nonzero);
    [X, genelist] = sc_rmdugenes(X, genelist);
    [X, keptidx] = sc_selectc(X, libszcutoff, gnnumcutoff);
    keptidxv{end+1} = keptidx;
    newsz = size(X);
    c = c + 1;
end
end


function [keptidx] = obj_rmmtcells(Xobj, genelist, mtratio, mtgenenamepat, vebrose)
%Remove cells with high mtDNA ratio
if nargin < 3, mtratio = 0.1; end
if nargin < 4, mtgenenamepat = "mt-"; end
if nargin < 5, vebrose = false; end

assert(isa(Xobj, 'pkg.refwrap'))
assert(size(Xobj.data, 1) == length(genelist))
idx = startsWith(genelist, mtgenenamepat, 'IgnoreCase', true);
if sum(idx) > 0
    if vebrose
        fprintf('%d mt-genes found.\n', sum(idx));
    end
    lbsz = sum(Xobj.data, 1);
    lbsz_mt = sum(Xobj.data(idx, :), 1);
    f_mtreads = lbsz_mt ./ lbsz;
    keptidx = f_mtreads < mtratio;
    if sum(~keptidx) > 0
        %Xobj.data=Xobj.data(:,keptidx);
        Xobj.data(:, ~keptidx) = [];
        if vebrose
            fprintf('%d cells with mt-read ratio >=%f (or %f%%) are removed.\n', ...
                sum(~keptidx), mtratio, mtratio*100);
        end
    else
        fprintf('No cells with mt-read ratio >=%f (or %f%%) are removed.\n', ...
            mtratio, mtratio*100);
    end
else
    if vebrose
        fprintf('No mt-genes found.\n');
        fprintf('No cells with mt-read ratio >=%f (or %f%%) are removed.\n', ...
            mtratio, mtratio*100);
    end
    if nargout > 1
        keptidx = true(size(Xobj.data, 2), 1);
    end
end
end

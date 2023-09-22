function [status] = sc_sce2hdf5(sce, filename)

% import h5py
% f = h5py.File("file.h5ad", "r")
% list(f.keys())
% f["sce"].visititems(print)
% f["sce/X"][()]


status = 0;
if nargin < 2
    [filen, pathn] = uiputfile( ...
        {'*.h5', '*.hdf5', '*.*'}, 'Save as');
    if ~(filename), return; end
    filename = [pathn, filen];
end

h5create(filename, '/sce/X', size(sce.X), 'datatype', 'uint16');
h5write(filename, '/sce/X', uint16(sce.X))
h5writeatt(filename, '/sce/X', 'source', "sc_sce2hdf5");

h5create(filename, '/sce/g', size(sce.g), 'Datatype', 'string')
h5write(filename, '/sce/g', sce.g)
h5create(filename, '/sce/s', size(sce.s))
h5write(filename, '/sce/s', sce.s)
h5create(filename, '/sce/c', size(sce.c))
h5write(filename, '/sce/c', sce.c)


% h5read(filename,'/matrix/data')

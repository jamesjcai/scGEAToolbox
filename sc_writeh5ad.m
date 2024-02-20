
data = rand(100, 10); % Example data, replace with your own
g="g"+string(1:100);
c="c"+string(1:10);

% py.scipy.sparse.csr_matrix(data)
% py.str(c(1))

filename = 'example_data.h5'; 
delete(filename)

% py.scipy.sparse.csr_matrix(data)

h5create(filename, '/X/data', size(data));
h5write(filename, '/X/data', py.numpy.array(data));

h5create(filename, '/obs/barcode', size(c(:)),'Datatype','string');
h5create(filename, '/var/gene', size(g(:)),'Datatype','string');

h5write(filename, '/obs/barcode', c(:));
h5write(filename, '/var/gene', g(:));

a = h5info("example_data.h5")



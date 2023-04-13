


% pkg.e_writeh5(X,genelist,filename)
h5create(filename, '/X', size(X));
h5write(filename, '/X', X);
h5create(filename, '/g', size(genelist),'Datatype','string');
h5write(filename, '/g', g);

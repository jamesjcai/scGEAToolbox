load ../../../example_data/new_example_sce.mat
X = sce.X;
if issparse(X)
    X = (full(X));
end
pyrun("import numpy as np")
pyrun("import scrublet as scr");

x = py.numpy.array(X);
pyrun("scrub = scr.Scrublet(X)",X=x);

pyrun("scrub = scr.Scrublet(np.array(X))",X=X);
a = pyrun("doubletscore, isDoublet = scrub.scrub_doublets()", ...
    "doubletscore");



pyrun("import h5py");
pyrun("import scrublet as scr");
pyrun("f = h5py.File('input.mat','r')");
pyrun("X = f['X'][()]");
pyrun("f.close()");
pyrun("scrub = scr.Scrublet(X)");
pyrun("doubletscore, isDoublet = scrub.scrub_doublets()");

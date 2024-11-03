function r_scplotter(sce, wkdir, isdebug)

if nargin < 3, isdebug = false; end
if nargin < 2, wkdir = tempdir; end

oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_scplotter');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5', 'output.png'};
if exist("input.h5", "file"), delete("input.h5"); end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%if ~strcmp(unique(sce.c_cell_type_tx), "undetermined")
    pkg.e_writeh5(full(sce.X), sce.g, 'input.h5', sce.c_cell_type_tx, sce.s);
%else
%    pkg.e_writeh5(full(sce.X), sce.g, 'input.h5');
%end

%sc_writefile('input.txt',sce.X,sce.g);
%    if isdebug, return; end
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);
if exist('output.png', 'file')
    img = imread('output.png');
    imshow(img);
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end


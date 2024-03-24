function [T] = r_SPIA(Tdeg, fdrcutoff, speciestag, wkdir)

if nargin < 4, wkdir = tempdir; end
if nargin < 2, fdrcutoff = 0.01; end
if nargin < 3, speciestag = 'hsa'; end
T = [];
if ~gui.gui_showrefinfo('SPIA [PMID:18990722]'), return; end

isdebug = false;
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_SPIA');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end


gid = pkg.i_symbol2ncbiid(Tdeg.gene);
Tdeg = Tdeg(gid ~= 0, :);
gid = gid(gid ~= 0);
writematrix(string(gid), 'input1.txt', 'QuoteStrings', 'all');

idx = Tdeg.p_val_adj < fdrcutoff;
gid = gid(idx);
Tdeg = Tdeg(idx, :);

writematrix(string(gid), 'input2.txt', 'QuoteStrings', 'all');
fc = Tdeg.avg_log2FC;
mfc = 2 * max(abs(fc(~isinf(fc))));
fc(isinf(fc)) = sign(fc(isinf(fc))) * mfc;
writematrix(fc, 'input3.txt');

tmpfilelist = {'input1.txt', 'input2.txt', 'input3.txt', 'output.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.RunRcode(codefullpath, Rpath);

if ~exist('output.csv', 'file'), return; end
T = readtable('output.csv');

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

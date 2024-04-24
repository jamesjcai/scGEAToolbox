function [T] = r_DESeq2(X, Y, genelist, wkdir)

if nargin < 4, wkdir = tempdir; end
if nargin < 3, genelist = (1:size(X, 1))'; end
T = [];
isdebug = false;
oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_DESeq2');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

if issparse(X), X = full(X); end
if issparse(Y), Y = full(Y); end

avg_1 = mean(X, 2);
avg_2 = mean(Y, 2);
pct_1 = sum(X > 0, 2) ./ size(X, 2);
pct_2 = sum(Y > 0, 2) ./ size(Y, 2);

%T = table(gene, p_val, avg_log2FC, abs_log2FC, avg_1, avg_2, ...
%    pct_1, pct_2, p_val_adj);

tmpfilelist = {'input.mat', 'output.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
if issparse(Y), Y = full(Y); end

save('input.mat', 'X', 'Y', '-v7.3');

Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepth,'script.R');
pkg.RunRcode(codefullpath, Rpath);

if ~exist('output.csv', 'file'), return; end
warning off
T = readtable('output.csv', 'TreatAsMissing', 'NA');
T.Var1 = genelist(T.Var1);
T.Properties.VariableNames{'Var1'} = 'gene';
T.Properties.VariableNames{'log2FoldChange'} = 'avg_log2FC';
abs_log2FC = abs(T.avg_log2FC);
T = addvars(T, abs_log2FC, 'After', 'avg_log2FC');

T = addvars(T, pct_2, 'After', 'abs_log2FC');
T = addvars(T, pct_1, 'After', 'abs_log2FC');
T = addvars(T, avg_2, 'After', 'abs_log2FC');
T = addvars(T, avg_1, 'After', 'abs_log2FC');

T = sortrows(T, 'abs_log2FC', 'descend');
T = sortrows(T, 'padj', 'ascend');
T.Properties.VariableNames{'padj'} = 'p_val_adj';
warning on
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

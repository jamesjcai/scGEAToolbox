function [T] = r_MAST(X, Y, genelist, wkdir)

if nargin < 4, wkdir = tempdir; end
if nargin < 3, genelist = (1:size(X, 1))'; end
T = [];
isdebug = false;
oldpth = pwd();
cleanupObj = onCleanup(@() cd(oldpth));

[isok, msg, codepath] = commoncheck_R('R_MAST');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

avg_1 = mean(X, 2);
avg_2 = mean(Y, 2);
% pct_1 = sum(X>0,2)./size(X,2);
% pct_2 = sum(Y>0,2)./size(Y,2);

tmpfilelist = {'input.mat', 'output.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
if issparse(X), X = full(X); end
if issparse(Y), Y = full(Y); end

save('input.mat', 'X', 'Y', '-v7.3');

Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.i_runrcode(codefullpath, Rpath);

if ~exist('output.csv', 'file'), return; end
T = readtable('output.csv', 'VariableNamingRule', 'modify');


% pct_2a=pct_2(T.Var1);
% pct_1a=pct_1(T.Var1);
avg_2 = avg_2(T.Var1);
avg_1 = avg_1(T.Var1);


T.Var1 = genelist(T.Var1);
T.Properties.VariableNames{'Var1'} = 'gene';
abs_log2FC = abs(T.avg_log2FC);
T = addvars(T, abs_log2FC, 'After', 'avg_log2FC');
%  T = addvars(T,pct_2a,'After','abs_log2FC');
%  T = addvars(T,pct_1a,'After','abs_log2FC');
T = addvars(T, avg_2, 'After', 'abs_log2FC');
T = addvars(T, avg_1, 'After', 'abs_log2FC');

T = sortrows(T, 'abs_log2FC', 'descend');
T = sortrows(T, 'p_val_adj', 'ascend');
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
% cd(oldpth);
end

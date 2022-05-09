function [T]=MAST(X,Y,genelist)

if nargin<3, genelist=(1:size(X,1))'; end
T=[];
isdebug=true;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_MAST');
if ~isok, error(msg); end

tmpfilelist={'input.mat','output.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
save('input.mat','X','Y','-v7.3');
pkg.RunRcode('script.R');
if ~exist('output.csv','file'), return; end
warning off
T=readtable('output.csv');
T.Var1=genelist(T.Var1);
T.Properties.VariableNames{'Var1'} = 'gene';
abs_log2FC=abs(T.avg_log2FC);
T = addvars(T,abs_log2FC,'After','avg_log2FC');
T=sortrows(T,'abs_log2FC','descend');
T=sortrows(T,'p_val_adj','ascend');    
warning on
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

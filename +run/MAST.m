function [T]=MAST(X,Y,genelist)
if nargin<3
    genelist=(1:size(X,1))';
end

oldpth=pwd();
[isok,msg]=commoncheck_R('R_MAST');
if ~isok, error(msg); end

if exist('output.csv','file')
    delete('output.csv');
end
writematrix(X,'input1.txt');
writematrix(Y,'input2.txt');
RunRcode('script.R');
if exist('output.csv','file')
    warning off
    T=readtable('output.csv');
    T.Var1=genelist(T.Var1);
    T.Properties.VariableNames{'Var1'} = 'gene';
    abs_log2FC=abs(T.avg_log2FC);
    T = addvars(T,abs_log2FC,'After','avg_log2FC');
    T=sortrows(T,'abs_log2FC','descend');
    T=sortrows(T,'p_val_adj','ascend');
    
    warning on
else
    T=[];
end
if exist('input1.txt','file'), delete('input1.txt'); end
if exist('input2.txt','file'), delete('input2.txt'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end

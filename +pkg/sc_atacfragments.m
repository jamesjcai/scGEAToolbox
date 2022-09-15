% function [X,g,b]=sc_atacfragments(infile)

infile='atac_fragments (1).tsv';
T=readtable(infile,'FileType','text','CommentStyle',"#", ...
                    'Format','%s%d%d%s%d');
[b,~,xid]=unique(string(T.Var4));
% peakid=strings(length(a.Var1),1);
% for k=1:length(peakid)
%     peakid(k)=sprintf('%s_%d_%d',a.Var1{k},a.Var2(k),a.Var3(k));
% end
% g=peakid;
g=strcat(T.Var1,"_",string(T.Var2),"_",string(T.Var3));
X=sparse((1:length(T.Var1))',xid,double(T.Var5));





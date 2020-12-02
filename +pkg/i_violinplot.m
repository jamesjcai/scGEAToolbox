function out=i_violinplot(d, c)
import pkg.Violin
import pkg.violinplot
if isstring(c)
    c=strrep(c,'_',' ');
end
[~,cL]=grp2idx(c);
[~,i]=sort(grpstats(d,c,@mean),'descend');
out=violinplot(d,c,'GroupOrder',cL(i),...
    'ShowData',false,'ViolinColor',[1 1 1],...
    'EdgeColor',[0 0 0],'ShowMean',false);
xtickangle(-45);
box on

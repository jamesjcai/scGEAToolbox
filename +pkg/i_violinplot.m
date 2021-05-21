function i_violinplot(d,c,colorit)
import pkg.Violin
import pkg.violinplot
% if isstring(c)
%     c=strrep(c,'_',' ');
% end
% [~,cL]=grp2idx(c);
% [~,i]=sort(grpstats(d,c,@median),'descend');
if nargin<3, colorit=false; end
if issparse(d), d=full(d); end
if ~colorit
    violinplot(d,c,...
        'ShowData',false,'ViolinColor',[1 1 1],...
        'EdgeColor',[0 0 0]);
else
    violinplot(d,c,'ShowData',false,'EdgeColor',[0 0 0]);
end
%xtickangle(-45);
box on
grid on
end

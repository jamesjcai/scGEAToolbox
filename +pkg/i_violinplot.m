function i_violinplot(d,c,colorit,grouporder)
import pkg.Violin
import pkg.violinplot
% if isstring(c)
%     c=strrep(c,'_',' ');
% end
% [~,cL]=grp2idx(c);
% [~,i]=sort(grpstats(d,c,@median),'descend');

if ~isstring(c)
    c=string(c);
end
if nargin<4, grouporder=[]; end
if nargin<3 || isempty(colorit), colorit=false; end
if issparse(d), d=full(d); end

if ~colorit
    if isempty(grouporder)
        violinplot(d,c,...
            'ShowData',false,'ViolinColor',[1 1 1],...
            'EdgeColor',[0 0 0]);
    else
        violinplot(d,c,...
            'ShowData',false,'ViolinColor',[1 1 1],...
            'EdgeColor',[0 0 0],'GroupOrder',grouporder);
    end
else
    if isempty(grouporder)
        violinplot(d,c,'ShowData',false,'EdgeColor',[0 0 0]);
    else
        violinplot(d,c,'ShowData',false,'EdgeColor',[0 0 0], ...
            'GroupOrder',grouporder);
    end
end
%xtickangle(-45);
box on
grid on
end

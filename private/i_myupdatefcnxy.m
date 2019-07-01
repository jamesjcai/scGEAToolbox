function txt = i_myupdatefcnxy(~,event_obj,g,F1,F2,X,Y)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
P1=F1(idx,:);
P2=F2(idx,:); 
pts = [P1; P2];
% line(pts(:,1), pts(:,2), pts(:,3))
plot3(pts(:,1), pts(:,2), pts(:,3),'k-')
return;

x=X(idx,:);
y=Y(idx,:);
persistent hFig
if ~isempty(hFig) && ishandle(hFig) && isvalid(hFig)
    figure(hFig)
else
    hFig = figure;
end
n=length(x);
m=length(y);
stem(1:n,x,'marker','none');
hold on
stem(n+1:n+m,y,'marker','none');
set(hFig,'Toolbar','none');
set(hFig, 'MenuBar', 'none');
set(hFig,'NumberTitle','off');
title(g(idx));
hold off
end
function txt = i_myupdatefcnx(~,event_obj,g,F1,F2)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
P1=F1(idx,:);
P2=F2(idx,:); 
pts = [P1; P2];
line(pts(:,1), pts(:,2), pts(:,3),'color','g');
legend({'Data 1','Spline 1','Data 2','Spline 2'},'location','northeast');
%plot3(pts(:,1), pts(:,2), pts(:,3),'k-')
% if ~isempty(findobj(gcf,'type','Legend'))
%     a=legend;
%     legend(a.String,'location',a.Location);
% end
end
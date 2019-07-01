function txt = i_myupdatefcn2(~,event_obj,g)
% Customizes text of data tips
pos = event_obj.Position;
% idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
% num2str(pos(2))
txt = {[num2str(pos(2)) ' - ' char(g(pos(2)))]};
% txt={num2str(pos(2))}
end
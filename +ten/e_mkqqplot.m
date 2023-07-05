function e_mkqqplot(T)

pd = makedist('Gamma','a',0.5,'b',2);
qqplot(T.FC,pd);

[~,idx]=sort(T.FC);
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1x,T.genelist(idx)};


% h1=h(1);
% h1.DataTipTemplate.DataTipRows = T.genelist(idx);
% for k=1:5
%     datatip(h1, 'DataIndex', idx(k));    
% end

end

function txt = i_myupdatefcn1x(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end

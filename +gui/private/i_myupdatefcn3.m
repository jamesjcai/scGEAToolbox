function txt = i_myupdatefcn3(~, event_obj, g, X, Y)
if nargin < 5
    Y = [];
end
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
persistent myupdatefcn3fig
if isempty(myupdatefcn3fig) || ~isvalid(myupdatefcn3fig)
    myupdatefcn3fig = figure;
    p = myupdatefcn3fig.Position;
    myupdatefcn3fig.Position = [p(1:3), 320];
    % myupdatefcn3fig.ToolBar='none';
    % myupdatefcn3fig.MenuBar='none';
end
if isvalid(myupdatefcn3fig) && isa(myupdatefcn3fig, 'matlab.ui.Figure')
    figure(myupdatefcn3fig);
end
stem(1:length(X(idx, :)), X(idx, :), 'marker', 'none');
if ~isempty(Y)
    hold on
    stem(1+length(X(idx, :)):length(Y(idx, :))+length(X(idx, :)), ...
        Y(idx, :), 'marker', 'none');
end
title(txt)
xlabel('Cell Index')
ylabel('Expression Level')
end    


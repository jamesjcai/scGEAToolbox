function txt = i_myupdatefcn4(~, event_obj, g, X, Y)
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
n1=length(X(idx, :));
x1 = X(idx, :);
stem(1:n1, x1, 'marker', 'none');
hold on
x2=-1*(x1==0);
stem(1:n1, x2, 'marker', 'none');

if ~isempty(Y)
    hold on
    stem(1+n1:length(Y(idx, :))+n1, ...
        Y(idx, :), 'marker', 'none');
end
title(txt)
xlabel('Cell Index')
ylabel('Expression Level')
end    


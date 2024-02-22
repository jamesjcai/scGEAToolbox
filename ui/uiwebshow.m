function uiwebshow(url)

if nargin < 1 || isempty(url), url = 'https://scgeatool.github.io/'; end


fig = uifigure();
h = uihtml(fig,'Position',[50 80 400 400]);
h.HTMLSource = webread(url);

% fig = uifigure('Name','title','Position',[100 100 500 500], ...
%     'CloseRequestFcn',@(fig,event) Cancel_Pushed(fig));

   

    % jOpenBrowser = javaObjectEDT(javax.swing.JButton('Open in browser'));
    % jOpenBrowser.setToolTipText('Open this webpage in system browser');
    % jhOpen = handle(jOpenBrowser, 'CallbackProperties');
    % set(jhOpen, 'ActionPerformedCallback', @(h,e)web(char(jAddressField.getText)));


OK = uibutton(fig,'Position',[50 30 100 22],'Text','OK','ButtonPushedFcn', @(OK,event) OK_Pushed(fig,table));
Cancel = uibutton(fig,'Position',[170 30 100 22],'Text','Cancel','ButtonPushedFcn', @(Cancel,event) Cancel_Pushed(fig));



function OK_Pushed(fig)
       delete(fig); % Close figure
end
function Cancel_Pushed(fig)
    delete(fig); % Close figure
end    

end
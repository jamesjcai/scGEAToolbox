function uiwebshow(url)

if nargin < 1 || isempty(url), url = 'https://scgeatool.github.io/'; end

g = ["INS1","APOE","H2-K1","PLP1","CRYAB","S100A6","H2-D1","GATM","DBI","VIM","FXYD1","LY6A","RARRES2","HBB-BS","LDHB"];

s = urlencode(strtrim(sprintf('%s\r',g)));
% https://string-db.org/api/image/network?identifiers=DRD1_HUMAN%0dDRD2_HUMAN&species=9606
url=sprintf('https://string-db.org/api/image/network?identifiers=%s&species=9606',s);
%web(url)

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
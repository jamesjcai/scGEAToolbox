function uitabgroupfig(hc,parentfig)

% f1=figure('visibl','off');
% a1=axes('Parent',f1);
% scatter(a1,randn(300,1),rand(300,1));
% f2=figure('visibl','off');
% a2=axes('Parent',f2);
% scatter(a2,randn(300,1),rand(300,1),'r');
% % figure; h2=scatter(randn(300,1),rand(300,1),'rs');
% uitabgroupfig({a1,a2})

if nargin<2, parentfig=[]; end
if ~isa(hc,'cell'), error('aaa'); end

n=length(hc);

hFigure=figure('Visible','off');
if ~isempty(parentfig)
    height = 420;
    width = 560;
    sz = parentfig.Position;
    x = sz(1);
    y = sz(2);
    %Position= [(x - width)/2, (y - height)/2, width, height];
    Position= [x, y, width, height];
    hFigure.Position=Position;
else
    movegui(hFigure,'center');
end
set(hFigure, 'MenuBar', 'none');
set(hFigure, 'ToolBar', 'none');

tabgp = uitabgroup(hFigure,"Position",[.0 .0 1.0 1.0]);

for k=1:n
    tag=sprintf('%d',k);
    tab1 = uitab(tabgp,"Title",tag,'Tag',tag);
    set(hc{k},'Parent',tab1);
    %hax1 = axes('Parent', tab1);
    %set(hc{k},'Parent',hax1);
    % tab2 = uitab(tabgp,"Title","Options",'Tag',"2");
    % tabgp.SelectionChangedFcn=@displaySelection;
end



%dt = datacursormode(hFigure);
%dt.UpdateFcn = {@i_myupdatefcnx};

drawnow;
set(hFigure,'Visible','on');


end


function displaySelection(src,event)
    t = event.NewValue;
    title = t.Title;
    disp("Viewing the " + title + " tab")
    % h1.Visible="on";
    % h2.Visible="off";
end




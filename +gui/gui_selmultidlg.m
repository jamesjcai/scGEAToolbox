function gui_selmultidlg(genelist)
txt_cell_array = {'line1';'line2';'line3';'line4';'line5'};

if nargin<1, genelist=txt_cell_array; end
if ~iscell(genelist)
    genelist=cellstr(genelist);
end    

f = figure;
% ax = axes(f);
% ax.Units = 'pixels';
% ax.Position = [75 75 325 280]

c = uicontrol('style','pushbutton','Position',[250 160 100 30],...
    'String','>','Callback', {@plotButtonPushed,genelist});

c = uicontrol('style','pushbutton',...
    'Position',[250 195 100 30]);
c.String = '<';
c.Callback = {@plotButtonPushed2,genelist};

h_list = uicontrol('style','list','max',length(genelist),...
   'min',1,'Position',[20 20 150 360],...
   'string',genelist);

h_list2 = uicontrol('style','list','max',length(genelist),...
   'min',1,'Position',[380 20 150 360],...
   'string',[]);

    function plotButtonPushed(src,event,genelist)
        if ~isempty(h_list.String)
        h_list2.String=...
            unique([h_list2.String; h_list.String(h_list.Value)]);        
        h_list.String=setxor(genelist,h_list2.String,'stable');
        
            set(h_list,'Value',1);
        end
    end

    function plotButtonPushed2(src,event,genelist)
        % https://www.mathworks.com/matlabcentral/answers/92064-why-do-i-receive-a-warning-when-i-repopulate-my-listbox-uicontrol-in-matlab
        % set(h_list2,'Value',1);
        if ~isempty(h_list2.String)
            h_list2.String(h_list2.Value)=[];        
            set(h_list2,'Value',1);
        end        
        h_list.String=setxor(genelist,h_list2.String,'stable');
    end
end
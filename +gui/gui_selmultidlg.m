function gui_selmultidlg(genelist)
% if ~iscell(genelist)
%     genelist=cellstr(genelist);
% end
txt_cell_array = {'line1';'line2';'line3';'line4';'line5'};
if nargin<1, genelist=txt_cell_array; end
    

f = figure;
% ax = axes(f);
% ax.Units = 'pixels';
% ax.Position = [75 75 325 280]
c = uicontrol('style','pushbutton','Position',[250 160 100 30]);
c.String = '>';
c.Callback = @plotButtonPushed;

c = uicontrol('style','pushbutton',...
    'Position',[250 195 100 30]);
c.String = '<';
c.Callback = @plotButtonPushed2;


h_list = uicontrol('style','list','max',20,...
   'min',1,'Position',[20 20 150 360],...
   'string',genelist);

h_list2 = uicontrol('style','list','max',20,...
   'min',1,'Position',[380 20 150 360],...
   'string',"");

    function plotButtonPushed(src,event)
        h_list2.String=...
            unique([h_list2.String; h_list.String(h_list.Value)]);
    end
    function plotButtonPushed2(src,event)
        % https://www.mathworks.com/matlabcentral/answers/92064-why-do-i-receive-a-warning-when-i-repopulate-my-listbox-uicontrol-in-matlab
        set(h_list2,'Value',1);
        if ~isempty(h_list2.String)
            h_list2.String(h_list2.Value)=[];
        end
    end

end
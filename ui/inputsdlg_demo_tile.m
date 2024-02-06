% demo script for inputsdlg_struct
%
% inputsdlg with struct as input allows to edit structure variables with size > 1

clear; close all;

Prompt = {'Name'
          'Vector'
          'Value'
          'String'
          'List'
          'Color'
          'Check'};

Formats = struct('type',{'edit','edit','edit','edit','list','color','check'});
Formats(1).enable = 'inactive'; % acts as a label
Formats(2).format = 'vector';
Formats(3).format = 'float';
Formats(5).style = 'popupmenu';
Formats(5).items = {'bluered','autumn','bone','colorcube','cool','copper','gray','hot','hsv','jet'};

Formats = Formats.';

% create initial answer
DefAns = cell(size(Prompt,1),5);
for datanr = 1:5
    % no editable text / FORMAT = 'text'
    DefAns{1,datanr} = ['Data ' num2str(datanr)];
    
    % numeric vector / FORMAT = 'vector' 
    DefAns{2,datanr} = randi(20,1,3);
    
    % numeric single value, no FORMAT needed
    DefAns{3,datanr} = rand(1);
    
    % editable text, no FORMAT needed
    DefAns{4,datanr} = char(randi([32 126],1,randi(20,1)));
    
    % popupmenu with strings / FORMAT = 'cellarray with strings'
    DefAns{5,datanr} = randi(numel(Formats(5).items),1);
    
    % color, no FORMAT needed
    DefAns{6,datanr} = rand(1,3);
    
    % logical value, no FORMAT needed
    DefAns{7,datanr} = logical(randi(1,1));
end

DefAns = cell2struct(DefAns,Prompt,1);
Prompt = repmat(Prompt,1,2);

Title = 'test';
Options.AlignControls = 'on';
Options.CreateFcn = @(~,~,handles)celldisp(get(handles,'type'));
Options.DeleteFcn = @(~,~,handles)celldisp(get(handles,'type'));

inputsdlg(Prompt,Title,Formats,DefAns,Options)
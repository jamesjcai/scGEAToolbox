%     prompt = {'Enter font size (1-30):'};
%     dlgtitle = 'Input';
%     dims = [1 35];
%     definput = {'15'};
%     if stxtyes, definput = {'10'}; end
%     answer = inputdlg(prompt,dlgtitle,dims,definput);
%     try
%         fsz=round(str2double(answer{1}));
%     catch
%         return;
%     end
%     if ~(fsz>=1 && fsz<=30), return; end

%         if stxtyes
%             stxt=sprintf('%s',cL{i});
%         else
%             stxt=sprintf('%d',i);
%         end
%         stxt=strrep(stxt,'_','\_');
%
%        if ~isempty(h.ZData)    % size(sce.s,2)==3
%            dtp=datatip(h,'DataIndex',idx(k)); % si(:,1),si(:,2),si(:,3),'SnapToDataVertex','off');
%             text(si(:,1),si(:,2),si(:,3),stxt,...
%                 'fontsize',fsz,'FontWeight','bold',...
%                 'BackgroundColor','w','EdgeColor','k');
%        else
%            dtp=datatip(h,si(:,1),si(:,2),'SnapToDataVertex','off');
%             text(si(:,1),si(:,2),stxt,...
%                 'fontsize',fsz,'FontWeight','bold',...
%                 'BackgroundColor','w','EdgeColor','k');
%        end

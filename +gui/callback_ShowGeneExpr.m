function callback_ShowGeneExpr(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [axx,bxx]=view();
    % if any([axx,bxx]==0), axx=ax; bxx=bx; end
    % gsorted=sort(sce.g);
    
answer = questdlg('Show expression of single or mulitple genes?',...
    'Single/Multiple Genes','Single','Multiple','Cancel','Single');

switch answer
    case 'Single'
        [gsorted]=gui.i_sortgenenames(sce);
        if isempty(gsorted), return; end
        [indx,tf] = listdlg('PromptString',{'Select a gene or select multiple genes to display individually','',''},...
            'SelectionMode','single','ListString',gsorted);
        if tf==1
            for k=1:length(indx)
                gui.i_cascadefig(sce,gsorted(indx(k)),axx,bxx,k);
            end
        end
    case 'Multiple'
        [glist]=gui.i_selectngenes(sce);
        if isempty(glist)
            helpdlg('No gene selected.','');
            return;
        %[gsorted]=gui.i_sortgenenames(sce);
        %if isempty(gsorted), return; end
        %[idx]=gui.i_selmultidlg(gsorted);
        %if isempty(idx), return; end
        %if isscalar(idx) && idx==0
        %   helpdlg('No gene selected.','');
        %    return;
        else
            [y,i]=ismember(upper(glist),upper(sce.g));
       if ~all(y), error('Unspecific running error.'); end
        %[~,i]=ismember(gsorted(idx),sce.g);
        x=sum(sce.X(i,:),1);
        if length(i)==1
           g=sce.g(i);
        elseif length(i)>1
            answer2=questdlg('Union (OR) or Intersection (AND)',...
                '','Individually',...
                'Intersection (AND)',...
                'Union (OR)',...                
                'Individually');
            switch answer2
                case 'Union (OR)'
                    g=sprintf("%s | ",glist); 
                case 'Intersection (AND)'
                    g=sprintf("%s & ",glist);
                    ix=sum(sce.X(i,:)>0,1)==length(i);
                    if ~any(ix)
                        helpdlg('No cells expressing all selected genes.','');
                        return;
                    end
                    x=x.*ix;
                case 'Individually'
                    for k=1:length(glist)
                        gui.i_cascadefig(sce,glist(k),axx,bxx,k);
                        % i_showcascade(sce,gsorted(idx(k)),axx,bxx,k);
                    end
                    return;
                otherwise
                    return;
            end
            g=extractBefore(g,strlength(g)-2);
        end
            f=figure('visible','off');
            [h1]=sc_scattermarker(x,g,sce.s,g,5);
            title(g);
            view(h1,axx,bxx);
            movegui(f,'center');
            set(f,'visible','on');                  
        end
    case 'Cancel'
        % helpdlg('Action cancelled.','');
        return;
end

end

% function i_showcascade(sce,g,axx,bxx,k)
%         f = figure('visible','off');
%         [h1]=sc_scattermarker(sce.X,sce.g,sce.s,g,5);
%         view(h1,axx,bxx);
%         % movegui(f,'center');        
%         P = get(f,'Position');
%         set(f,'Position',[P(1)-20*k P(2)-20*k P(3) P(4)]);
%         set(f,'visible','on');
% end                       



% function [gsorted]=i_sortg(sce)
%         gsorted=[];
%         answer2 = questdlg('How to sort gene names?','Sort by',...
%             'Alphabetic','Expression Mean','Dropoff Rate','Alphabetic');
%         switch answer2
%             case 'Alphabetic'
%                 gsorted=sort(sce.g);
%             case 'Expression Mean'
%                 [T]=sc_genestats(sce.X,sce.g);
%                 [~,idx]=sort(T.Dropout_rate);
%                 gsorted=sce.g(idx);                
%             case 'Dropoff Rate'
%                 [T]=sc_genestats(sce.X,sce.g);
%                 [~,idx]=sort(T.Dropout_rate);
%                 gsorted=sce.g(idx);
%             otherwise
%                 return;
%         end
% end

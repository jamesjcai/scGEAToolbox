function callback_ShowGeneExpr(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [axx,bxx]=view();
    % if any([axx,bxx]==0), axx=ax; bxx=bx; end
    gsorted=sort(sce.g);
    
answer = questdlg('Show expression of single or mulitple genes?',...
	'Single/Multiple Genes',...
	'Single','Multiple','Cancel','Single');

switch answer
    case 'Single'
        [indx,tf] = listdlg('PromptString',{'Select a gene',...
        '',''},'SelectionMode','multiple','ListString',gsorted);
        if tf==1
            for k=1:length(indx)
                i_show(sce,gsorted(indx(k)),axx,bxx);
            end
        end
    case 'Multiple'
        [idx]=gui.gui_selmultidlg(gsorted);
        if isempty(idx), return; end
        if isscalar(idx) && idx==0
            helpdlg('No gene selected.');
            return;
        else
        [~,i]=ismember(gsorted(idx),sce.g);
        x=sum(sce.X(i,:),1);
        if length(i)==1
           g=sce.g(i);
        elseif length(i)>1
            answer2=questdlg('Union (OR) or Intersection (AND)','','Union (OR)','Intersection (AND)','Intersection (AND)');
            switch answer2
                case 'Union (OR)'
                    g=sprintf("%s | ",gsorted(idx)); 
                case 'Intersection (AND)'
                    g=sprintf("%s & ",gsorted(idx));                
                    ix=sum(sce.X(i,:)>0,1)==length(i);
                    if ~any(ix)
                        helpdlg('No cells expressing all selected genes.');
                        return;
                    end
                    x=x.*ix;
                otherwise
                    return;
            end            
            g=extractBefore(g,strlength(g)-2);
        end
            f = figure('visible','off');
            [h1]=sc_scattermarker(x,g,sce.s,g,5);
            title(g);
            view(h1,axx,bxx);
            movegui(f,'center');
            set(f,'visible','on');                  
        end
    case 'Cancel'
        % helpdlg('Action cancelled.');
        return;
end

end

function i_show(sce,g,axx,bxx)
        f = figure('visible','off');
        [h1]=sc_scattermarker(sce.X,sce.g,...
               sce.s,g,5);
        view(h1,axx,bxx);
        movegui(f,'center');
        set(f,'visible','on');
end

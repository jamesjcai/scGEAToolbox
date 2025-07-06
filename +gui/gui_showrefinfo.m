function [y, txt, T] = gui_showrefinfo(reftarget, parentfig)
    %see also: gui.gui_uishowrefinfo
    if nargin<2, parentfig = []; end
    
    y=false;
    txt = [];
    pw1 = fileparts(mfilename('fullpath'));
    fname = fullfile(pw1, '..','assets','Misc','refinfo.txt');
    fid=fopen(fname,'r');
    T=textscan(fid,'%s%s','Delimiter','\t');
    fclose(fid);
    reftag=string(T{:,1});
    
    idx=find(reftarget==reftag);
    if ~isempty(idx)
        txt=T{:,2}{idx};
    end
    
    if isempty(txt) , return; end
    fprintf('%s\n%s\n', reftarget, txt);
    
    if gui.i_isuifig(parentfig)
    
        hFig = uifigure("WindowStyle","modal",'Visible','off');
        hFig.Position(3)=0.75*hFig.Position(3);
        hFig.Position(4)=0.75*hFig.Position(4);
        g = uigridlayout(hFig,[3 3]);
        g.RowHeight = {'fit','2x','fit'};
        g.ColumnWidth = {'1x',75,75};
        
        lbl = uilabel(g,"Text",reftarget);
        lbl.Layout.Row = 1;
        lbl.Layout.Column = 1;
        
        txa = uitextarea(g);
        txa.Layout.Row = 2;
        txa.Layout.Column = [1 3];
        txa.Value = compose(txt);
        
        % txa.Position(4)=txa.Position(4)*2;
        if nargout>0
            oktext = "Continue";
        else
            oktext = "OK";
        end
        
        btn = uibutton(g,"Text",oktext);
        btn.Layout.Row = 3;
        btn.Layout.Column = 2 + (nargout==0);
        btn.ButtonPushedFcn = @(src,event) textEntered(src,event,btn);
        %btn.Position(3) = 50;
        
        if nargout > 0
            btn2 = uibutton(g,"Text","Cancel");
            btn2.Layout.Row = 3;
            btn2.Layout.Column = 3;
            btn2.ButtonPushedFcn = @(src,event) textEntered(src,event,btn2);
        end
        
        gui.i_movegui2parent(hFig, parentfig);
        
        % hFig.WindowStyle="modal";
        % pause(0.5);
        drawnow;
        hFig.Visible=true;
        % uiwait(hFig);
    
    else 
        answer = gui.myQuestdlg(parentfig, ...
             txt, reftarget, {'Continue','Cancel'}, 'Continue');
         switch answer
             case 'Continue'
                 y = true;
             otherwise
         end    
    end
    
    function textEntered(~,~,btn)
        if strcmp(btn.Text,oktext)
            y = true;
        end
        delete(hFig);
    end
    
end


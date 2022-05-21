function i_scintimeres(s,g)

FigureHandle=figure();
hAx = axes('Parent', FigureHandle);
c=ones(size(s,1),1);
[h] = gui.i_gscatter3(s, c, 1,1,hAx);

kc = numel(unique(c));
colormap(pkg.i_mycolorlines(kc));
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcnx};

tb = uitoolbar('Parent', FigureHandle);
pkg.i_addbutton2fig(tb,0,@refreshc,'aa.gif','clustering');
pkg.i_addbutton2fig(tb,0,@highlightmax,'aa.gif','clustering');
pkg.i_addbutton2fig(tb,0,@saveout,'aa.gif','clustering');

    function saveout(~,~)
%         txt="";
%         for k=1:max(c)
%             txt=sprintf('%s\t%s', ...
%                 txt,sprintf('%s ',g(c==k)));
%         end
        tmpName=[tempname,'.txt'];
        fid=fopen(tmpName,'w');
        for k=1:max(c)
            a=sprintf('%s ',g(c==k));
            fprintf(fid,'Cluster %d\t%s\n',k, a);
        end
        %fprintf(fid,'%s',txt);
        fclose(fid);
        [status]=system(['notepad "' tmpName '" &']);
        if status~=0
           edit(tmpName);
        end
    end


    function refreshc(~,~)
        kc = gui.i_inputnumk(50);
        if isempty(kc), return; end
        c=sc_cluster_s(s,kc,'type','kmeans');
        h = gui.i_gscatter3(s, c, 1, hAx);
    end

    function highlightmax(~,~)
        a=grpstats(s,c,@mean);
        d=sum(a.^2,2);
        [~,idx]=sort(d,'descend');
        c2=ones(size(c));
        for k=1:10
            c2(c==idx(k))=2;
        end
        h = gui.i_gscatter3(s, c2, 1, hAx);
    end

    function [txt] = i_myupdatefcnx(~, event_obj)
        % pos = event_obj.Position;
        idx = event_obj.DataIndex;
        txt = g(idx);        
    end

end


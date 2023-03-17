function [f]=i_violinplot(y,thisc,ttxt,colorit,cL)
    if nargin<5
        [~,cL]=grp2idx(thisc);
    end
    if nargin<4
        colorit=true;
    end
    f=figure('visible','off');
    tb=uitoolbar(f);
    pkg.i_addbutton2fig(tb,'off',{@i_savedata,y,thisc}, ...
        'export.gif','Export data...');
    pkg.i_addbutton2fig(tb,'off',{@i_testdata,y,thisc}, ...
        'exportx.gif','ANOVA/T-test...');    
    pkg.i_addbutton2fig(tb,'off',@i_addsamplesize, ...
        "xpowerpoint.gif",'Add Sample Size');    
    pkg.i_addbutton2fig(tb,'on',{@gui.i_savemainfig,3}, ...
        "powerpoint.gif",'Save Figure to PowerPoint File...');
    pkg.i_addbutton2fig(tb,'off',@i_invertcolor, ...
        "xpowerpoint.gif",'Switch BW/Color');
    pkg.i_addbutton2fig(tb,'off',@i_reordersamples, ...
        "xpowerpoint.gif",'Reorder Samples');  
    pkg.i_addbutton2fig(tb,'off',@i_sortbymean, ...
        "xpowerpoint.gif",'Sort Samples by Median');  

    isdescend=false;
    OldTitle=[];
    OldXTickLabel=[];
    cL=strrep(cL,'_','\_');
    thisc=strrep(string(thisc),'_','\_');
    pkg.i_violinplot(y,thisc,colorit,cL);
    title(strrep(ttxt,'_','\_'));
    %ylabel(selitems{indx1});
    movegui(f,'center');
    set(f,'visible','on');
    
%catch ME
%    errordlg(ME.message);
%end

    function i_invertcolor(~,~)
        colorit=~colorit;
        cla;
        pkg.i_violinplot(y,thisc,colorit,cL);
    end

    function i_addsamplesize(~,~)
        b=gca;
        if isempty(OldXTickLabel)
            a=zeros(length(cL),1);
            
            OldXTickLabel=b.XTickLabel;
            for k=1:length(cL)
                a(k)=sum(thisc==cL(k));
                b.XTickLabel{k}=sprintf('%s\\newline(n=%d)', ...
                    b.XTickLabel{k},a(k));
            end
        else
            b.XTickLabel=OldXTickLabel;
            OldXTickLabel=[];
        end
    end

    function i_sortbymean(~,~)
        [cx,cLx]=grp2idx(thisc);
        a=zeros(max(cx),1);
        for k=1:max(cx)
            a(k)=median(y(cx==k));
        end
        if isdescend
            [~,idx]=sort(a,'ascend');
            isdescend=false;
        else
            [~,idx]=sort(a,'descend');
            isdescend=true;
        end
        cLx_sorted=cLx(idx);

        %[~,cL,noanswer]=gui.i_reordergroups(thisc,cLx_sorted);
        %if noanswer, return; end
        cla
        pkg.i_violinplot(y,thisc,colorit,cLx_sorted);
    end


    function i_reordersamples(~,~)

        [~,cL,noanswer]=gui.i_reordergroups(thisc);
        if noanswer, return; end
        cla
        pkg.i_violinplot(y,thisc,colorit,cL);
    end


function i_testdata(~,~,y,grp)
    a=gca;
    if isempty(OldTitle)
        OldTitle=a.Title.String;
        if size(y,2)~=length(grp)
            y=y.';
        end
        tbl=pkg.e_grptest(y,grp);
    %h1=gca;
    %titre=string(h1.Title.String);
    
%     a=sprintf('%s\n%s=%.2e; %s=%.2e', ...
%         strrep(string(ttxt),'_','\_'), ...
%         strrep(tbl.Properties.VariableNames{1},'_','\_'), ...
%         tbl.(tbl.Properties.VariableNames{1}), ... 
%         strrep(tbl.Properties.VariableNames{2},'_','\_'), ...
%         tbl.(tbl.Properties.VariableNames{2}));
%     title(a);

    b=sprintf('%s=%.2e; %s=%.2e', ...
        strrep(tbl.Properties.VariableNames{1},'_','\_'), ...
        tbl.(tbl.Properties.VariableNames{1}), ... 
        strrep(tbl.Properties.VariableNames{2},'_','\_'), ...
        tbl.(tbl.Properties.VariableNames{2}));
        if iscell(OldTitle)
            newtitle=OldTitle;
        else
            newtitle={OldTitle};
        end
        newtitle{2}=b;
        a.Title.String=newtitle;
    else
        a.Title.String=OldTitle;
        OldTitle=[];
    %gui.i_exporttable(tbl,true);
    end
end


end

function i_savedata(~,~,a,b)
    T=table(a(:),b(:));    
    T.Properties.VariableNames={'ScoreLevel','GroupID'};
    T=sortrows(T,'ScoreLevel','descend');
    T=sortrows(T,'GroupID');
    gui.i_exporttable(T,true);
end



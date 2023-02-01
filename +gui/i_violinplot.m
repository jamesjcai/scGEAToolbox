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
    pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3}, ...
        "powerpoint.gif",'Save Figure to PowerPoint File...');
    pkg.i_addbutton2fig(tb,'off',@i_invertcolor, ...
        "xpowerpoint.gif",'Switch BW/Color');
    pkg.i_addbutton2fig(tb,'off',@i_reordersamples, ...
        "xpowerpoint.gif",'Reorder Samples');   

    cL=strrep(cL,'_','\_');
    thisc=strrep(thisc,'_','\_');
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


    function i_reordersamples(~,~)
        [~,cL,noanswer]=gui.i_reordergroups(thisc);
        if noanswer, return; end
        cla
        pkg.i_violinplot(y,thisc,colorit,cL);
    end

end

function i_savedata(~,~,a,b)
    T=table(a(:),b(:));    
    T.Properties.VariableNames={'ScoreLevel','GroupID'};
    T=sortrows(T,'ExprLevel','descend');
    T=sortrows(T,'GroupID');
    gui.i_exporttable(T,true);
end

function i_testdata(~,~,y,grp)
    if size(y,2)~=length(grp)
        y=y.';
    end
    tbl=pkg.e_grptest(y,grp);
    gui.i_exporttable(tbl,true);
end


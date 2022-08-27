function i_cascadeviolin(sce,Xt,thisc,glist,ytxt,grouporder,colorit)
if nargin<7, colorit=false; end
if nargin<6, grouporder=[]; end
if nargin<5, ytxt=''; end

F=cell(length(glist),1);
for k=1:length(glist)
    [~,idx]=ismember(glist(k),sce.g);
    y=full(Xt(idx,:));
    ttxt=sce.g(idx);
    
    f = figure('visible','off');
    pkg.i_violinplot(y,thisc,colorit,grouporder);
    title(strrep(ttxt,'_','\_'));
    ylabel(ytxt);
    tb=uitoolbar(f);

    pkg.i_addbutton2fig(tb,'off',{@i_savedata,y,thisc}, ...
        'export.gif','Export data...');
    pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3}, ...
        "powerpoint.gif",'Save Figure to PowerPoint File...');

    P = get(f,'Position');
    set(f,'Position',[P(1)-20*k P(2)-20*k P(3) P(4)]);
    set(f,'visible','on');                
    drawnow;
    F{k}=f;
end
gui.i_export2pptx(F,glist);
end

function i_savedata(~,~,a,b)
    T=table(a(:),b(:));    
    T.Properties.VariableNames={'ExprLevel','GroupID'};
    T=sortrows(T,'ExprLevel','descend');
    T=sortrows(T,'GroupID');
    gui.i_exporttable(T,true);
end

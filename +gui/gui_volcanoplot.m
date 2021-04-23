function gui_volcanoplot(T)
x=T.avg_logFC;
y=-log10(T.p_val_adj);
y(y>100)=100;
isok=T.pct_1>0.01|T.pct_2>0.01;
genelist=T.gene;

% figure;
ix=isinf(x) | abs(x)>10;
capv=1.05*max(abs(x(~ix)));
isx=sign(x);
x(ix)=capv;
x=isx.*abs(x);

scatter(x(isok),y(isok))
hold on
scatter(x(~isok),y(~isok),'x')
ylabel('-log10(adj p-value)')
xlabel('log2(FC)')
yline(-log10(0.01),'r')
xline(-1,'r')
xline(1,'r')

        dt=datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist};
        
        
function txt = i_myupdatefcn1(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end

end
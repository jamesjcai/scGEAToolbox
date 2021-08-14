function i_gjoyplot(X,genelist,targetg,g)

[c,cL]=grp2idx(g);
i=targetg==genelist;
x=X(i,:);
D=zeros(max(c),100);
for k=1:max(c)
    D(k,:)=ksdensity(X(i,c==k));
end
gui.i_joyplot(D,0.3,cL);


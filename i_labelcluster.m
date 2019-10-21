function i_labelcluster(s,c)
hold on
id=grp2idx(c);
for k=1:max(id)
    i=id==k;
    txtc=sprintf('%s (%d)',...
        unique(c(i)),numel(c(i)));    
    posc=mean(s(i,:));
    if length(posc)==2
        text(posc(1),posc(2),txtc);
    else
        text(posc(1),posc(2),posc(3),txtc);
    end       
end

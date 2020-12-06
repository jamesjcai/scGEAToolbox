function callback_PickColormap(~,~,kc)
    cx=colormap('autumn');
    cx(1,:)=[.8 .8 .8];
    a=lines(kc);
    b=a(randperm(size(a,1)),:);
    co={cx,a,b,'default','summer','jet','copper','winter'};
    colormap(co{randi(length(co))});
end
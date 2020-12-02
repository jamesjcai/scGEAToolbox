function PickColormap(~,~)
    cx=colormap('autumn');
    cx(1,:)=[.8 .8 .8];
    co={cx,'default','summer','jet','copper','winter'};
    colormap(co{randi(length(co))});
end
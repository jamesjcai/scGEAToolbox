function callback_PickColormap(~,~)
    cx=colormap('autumn');
    cx(1,:)=[.8 .8 .8];
    %a=lines(kc);
    %b=a(randperm(size(a,1)),:);
    mycmap=pkg.i_mycolormap;
    co={cx,'lines','default','summer',...
        'jet','copper','winter','hsv',mycmap};
%     co={cx,a,b,'default',summer(kc),...
%         jet(kc),copper(kc),winter(kc)};
    colormap(co{randi(length(co))});
end
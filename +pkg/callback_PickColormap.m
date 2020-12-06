function callback_PickColormap(~,~)
    cx=colormap('autumn');
    cx(1,:)=[.8 .8 .8];
    %a=lines(kc);
    %b=a(randperm(size(a,1)),:);
    co={cx,'autumn','lines','default','summer',...
        'jet','copper','winter'};
    
%     co={cx,a,b,'default',summer(kc),...
%         jet(kc),copper(kc),winter(kc)};
    colormap(co{randi(length(co))});
end
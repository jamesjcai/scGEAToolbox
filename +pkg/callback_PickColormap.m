function callback_PickColormap(~,~,pw1,n)
    if nargin<4, n=21; end
    if nargin<3, pw1=[]; end
    % disp(sprintf('Using %d colors',n));
    if n<3, n=3; end
    if ~isempty(pw1)
        pth=fullfile(pw1,'thirdparty/cbrewer');
        addpath(pth);
        CT=cbrewer('seq','Blues',n);
    end

    cx=autumn(n);
    cx(1,:)=[.8 .8 .8];
    %a=lines(kc);
    %b=a(randperm(size(a,1)),:);
    mycmap=pkg.i_mycolormap(n);
    co={cx,lines(n),'default',summer(n),...
        jet(n),copper(n),winter(n),hsv(n),...
        mycmap,CT,...
        cbrewer('div','Spectral',n),...
        cbrewer('div','RdBu',n),...
        cbrewer('seq','PuBuGn',n),...
        cbrewer('qual','Set1',n),...
        cbrewer('qual','Dark2',n)};
%     co={cx,a,b,'default',summer(kc),...
%         jet(kc),copper(kc),winter(kc)};
    colormap(co{randi(length(co))});
end
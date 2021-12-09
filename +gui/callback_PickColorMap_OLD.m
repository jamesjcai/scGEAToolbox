function callback_PickColorMap_OLD(~,~,n,showzero)
    if nargin<4, showzero=false; end
    if nargin<3, n=21; end
    % disp(sprintf('Using %d colors',n));
    n=max([n 3]);
    folder=fileparts(mfilename('fullpath'));
    a=strfind(folder,filesep);
    folder=extractBefore(folder,a(end)+1);
    wrkpth=fullfile(folder,'+run','thirdparty','cbrewer');
    addpath(wrkpth);
    CT=cbrewer('seq','Blues',n);

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
    if showzero
        cm=colormap;
        cm(1, :) = [.8 .8 .8];
        colormap(cm);
    end
end
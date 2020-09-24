function i_joyplot(D,scalek)

% James Cai (jcai@tamu.edu)
% (c) Aug 2017

    if nargin<2 || isempty(scalek), scalek=0.3; end
    D=D./max(D(:));
    [n,p]=size(D);
    hold on
    for k=n:-1:1
        %plot(1:p,scalek*(k-1)*ones(1,p)+D(k,:),'k-');
        i_joy2(D,k,p,scalek)
    end    
    set(gca, 'YTick', []);
    ylim([0 (n*scalek)+(1-scalek)]);
    xlim([1 p]);
    axis off
end

function i_joy2(D,k,p,scalek)
       y=scalek*(k-1)*ones(1,p)+D(k,:);
       y(1)=scalek*(k-1);
       y(end)=scalek*(k-1);
       x=1:p;
       patch(x,y,'red','FaceAlpha',.5);
end

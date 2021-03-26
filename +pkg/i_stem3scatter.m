function i_stem3scatter(x,y,c,t)
if nargin<4, t=''; end
figure;
    stem3(x,y,c,'marker','none','color','m');
    hold on
    scatter3(x,y,zeros(size(y)),5,c,'filled');
    title(t);
end
    

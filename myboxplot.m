function myboxplot(y,g,colorid)
if nargin<3
    colorid=1;
end
    c1=[0.70196 0.78039 1]; 
    c2=[1 0.6 0.78431]; 
    c3=[.8 .8 .8];
    c2=defaultcolor(2);
 switch colorid
     case 1
         c=c1;
     case 2
         c=c2;
     case 3
         c=c3;
     otherwise
         c='k';
 end
        %defaultc=get(gca,'ColorOrder');
        defaultc=get(groot,'DefaultAxesColorOrder');
        c=defaultc(colorid,:);
        %plot(1+g+0.1*(rand(size(g))-0.5),y,'o','linesmooth','off','color',c)
        %plot(1+g+0.1*(rand(size(g))-0.5),y,'o','color',c)
        
        %hold on
        g=double(g);
        h2=boxplot(y,g,'colors','k');
        set(h2(6,:),'color','k','linewidth',3);
        
        hold on
        % plot(g+1+randn(size(g))*0.025,d2,'o','color','k');
        %pp1=plot(1+g+0.1*(rand(size(g))-0.5),y,'o','color',c);
        %pp1=scatter(1+g+0.1*(rand(size(g))-0.5),y);
        [g]=i_add_jitter(g,y);
        %pp1=scatter(1+g+0.025*(randn(size(g))),y);
        pp1=scatter(g,y);
          %pp1.MarkerFaceColor=c;
        pp1.MarkerEdgeColor='k';
          %pp1.MarkerFaceAlpha=0.3;
        
        %set(h2,'linesmooth','off','linewidth',1.5);

        box on
        hold off
        
        % http://stackoverflow.com/questions/21999451/how-to-get-the-values-of-the-outliers-and-their-coordinates-from-a-box-plot
function [g]=i_add_jitter(g,A)
    gx=unique(g);
    for k=1:length(gx)
       i=g==gx(k);
       g(i)=i_add_j(g(i),A(i));
    end
    end

end
function [g]=i_add_j(g,A)        
    Q1 = quantile(A,0.25);
    Q3 = quantile(A,0.75);
    Spread = 1.5*(Q3-Q1);
    MaxValue = Q3 + Spread;
    MinValue = Q1 - Spread;
    i=~( A>MaxValue | A<MinValue);
    g(i)=g(i)+0.025*(randn(size(g(i))));
end

function c=defaultcolor(id)
% if nargin<1
%     c=get(gca,'ColorOrder');
%     return;
% end
%     c=get(gca,'ColorOrder');
%     c=c(id,:);
c=[      0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
if nargin==1
    c=c(id,:);
end
end


%function i_dotplot(X0,X1,genelist,tgene,uselog)
function i_dotplot(X,g,c,cL,tgene,uselog)

if nargin<6, uselog=false; end
[yes]=ismember(tgene,g);
if ~any(yes), return; end
z=length(tgene)-sum(yes);
if z>0
    fprintf('%d gene(s) not in the list are excluded.\n',z); 
end
tgene=tgene(yes);

% tgene=string(T.gene(1:10));
%idx=(1:length(tgene))';
%x=[-ones(size(idx)); ones(size(idx))]./2;
%y=repmat(idx,length(cL),1);

l=ones(length(tgene)*length(cL),1);
sz=l; vl=l;
x=l; y=l;
ct=0;
for k=1:length(tgene)
    for kk=1:length(cL)
        ct=ct+1;
        x(ct)=kk; y(ct)=k;
        a0=X(g==tgene(k),c==kk);
        sz(ct)=sum(a0~=0)./length(a0);
        vl(ct)=mean(a0);
    end
end


if uselog
    vl=log2(vl+1);
end
txgene=[" "; tgene];

% figure;
%sz=randi(100,1,length(x));
%scatter([-.5 .5],[-1 -1],[1 500],'k','filled');
%hold on
sz=sz+0.001;
vl=vl+0.001;
scatter(x,y,500*sz,vl,'filled');
hold on
scatter(x,y,500*sz,'k');
xlim([0.5 length(cL)+0.5]);
ylim([0.5 length(txgene)-0.5]);
colorbar
set(gca,'YTick',0:length(tgene))
set(gca,'YTickLabel',txgene)
set(gca,'XTick',0:length(cL))
set(gca,'XTickLabel',[{''};cL(:);{''}])
colormap(flipud(bone));
box on
grid on
hFig=gcf;
hFig.Position(3)=hFig.Position(3)*0.7;

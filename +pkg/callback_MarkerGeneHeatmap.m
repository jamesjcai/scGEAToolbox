function callback_MarkerGeneHeatmap(src,~)
sce=guidata(src);
[c,cL]=grp2idx(sce.c_cell_type_tx);
X=[]; sz=[]; cc=[];
for k=1:length(cL)
    i=c==k;
    X=[X sce.X(:,i)];
    sz=[sz sum(i)];
    cc=[cc; c(i)];
end

M=cell(numel(cL),2);
for k=1:numel(cL)
    [markerlist]=sc_pickmarkers(sce.X,sce.g,c,k);
    cLk=matlab.lang.makeValidName(cL{k});
    M{k,1}=cLk;
    M{k,2}=markerlist(1:50);    
end

% =========== 

subi=1:10:size(X,2);
for k=1:numel(cL)
    markerlist=M{k,2};
    cLk=M{k,1};
    [~,idx_g]=ismember(upper(markerlist),upper(sce.g));
    xx=X(idx_g,subi);
    ccc=cc(subi);
    
    M=[m_ncsc; m_csc; m_tac; m_ec; m_eec; m_goblet; m_dcs; m_tuft];
% M=unique(M,'stable');
[~,idx_g]=ismember(upper(M),upper(g));
X=Xori(idx_g,:);

X=sc_norm(X);
X=log2(X+1);


y=[];
c=[];
sz=[length(m_ncsc) length(m_csc) length(m_tac) length(m_ec) length(m_eec) length(m_goblet) length(m_dcs) length(m_tuft)];
i=sce.c_cell_type_tx=="Noncycling SC";
y=[y X(:,i)]; c=[c ones(1,sum(i))];
i=sce.c_cell_type_tx=="Cycling SC";
y=[y X(:,i)]; c=[c 2*ones(1,sum(i))];
i=sce.c_cell_type_tx=="TA";
y=[y X(:,i)]; c=[c 3*ones(1,sum(i))];
i=sce.c_cell_type_tx=="EC";
y=[y X(:,i)]; c=[c 4*ones(1,sum(i))];
i=sce.c_cell_type_tx=="EEC";
y=[y X(:,i)]; c=[c 5*ones(1,sum(i))];
i=sce.c_cell_type_tx=="Goblet (type 1)" | sce.c_cell_type_tx=="Goblet (type 2)";
y=[y X(:,i)]; c=[c 6*ones(1,sum(i))];
i=sce.c_cell_type_tx=="DCS (type 1)" | sce.c_cell_type_tx=="DCS (type 2)";
y=[y X(:,i)]; c=[c 7*ones(1,sum(i))];
i=sce.c_cell_type_tx=="Tuft cell";
y=[y X(:,i)]; c=[c 8*ones(1,sum(i))];



i=1:20:size(y,2);
cc=c(i);
xx=y(:,i);

figure;
% xx=quantilenorm(xx')';
    xx=zscore(xx,0,2);
    qx=quantile(xx(:),0.90);
    xx(xx>qx)=qx;
    qx=quantile(xx(:),0.10);
    xx(xx<qx)=qx;    
    imagesc(xx);

% figure;
% pkg.heatmap(xx, M,...
%     'ShowAllTicks', true);
% set(gca,'YTick',1:size(xx,1));

% axis xy
szc=cumsum(sz);
for k=1:max(cc)-1
    xline(sum(cc<k+1)+0.5,'r-');
    yline(szc(k)+0.5,'r-');
end
% a=colormap('autumn');
% a(1,:)=[.8 .8 .8];
% colormap(a);
set(gca,'YTick',1:size(xx,1));

set(gca,'XTick',[])
% set(gca,'YTick',[])
set(gca,'YTickLabel',M);
set(gca,'TickLength',[0 0])


%{    
    
    FigureHandle=figure;
    hAx = axes('Parent',FigureHandle);
    UitoolbarHandle = uitoolbar( 'Parent', FigureHandle ) ; 
    pt3 = uipushtool(UitoolbarHandle,'Separator','off');
    pt3.Tooltip = 'Select a gene to show expression';
    
    cLk=matlab.lang.makeValidName(cL{k});
    pt3.ClickedCallback = {@callback_PickMarkerGene,...
        markerlist,cLk};
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','greencircleicon.gif'));
ptImage = ind2rgb(img,map);

pt3.CData = ptImage;
    
    xxx=zscore(xx,0,2);
    qx=quantile(xxx(:),0.90);
    xxx(xxx>qx)=qx;

    %qx=quantile(xx(:),0.95);
    %xx(xx>qx)=qx;
%{    
xxk=xx;
try
for ix=1:size(xx,2)
    thisc=xx(ix,:);    
    [~,thisidx]=sort(thisc);
    thisa=ksdensity(thisc,"NumPoints",size(xx,2));
    thisc(thisidx)=thisa;
    xxk(ix,:)=thisc;
end
catch ME
    disp(ME)
end
%}    
    imagesc(xxx)
    
    for kk=1:max(ccc)-1
        xline(sum(ccc<kk+1)+0.5,'w-');
        % yline(szc(kk)+0.5,'w-');
    end
    set(gca,'YTick',1:size(xx,1));
    set(gca,'YTickLabel',markerlist);
    title(strrep(cLk,'_','\_'));
end

%}

end
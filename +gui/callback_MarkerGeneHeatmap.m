function callback_MarkerGeneHeatmap(src,~)

    answer = questdlg('Generate marker gene heatmap',...
        'Select Method','Method 1 (fast)','Method 2 (slower)',...
        'Method 3 (slowest)','Method 1 (fast)');
    switch answer
        case 'Method 1 (fast)'
            methodid=1;
        case 'Method 2 (slower)'
            methodid=3;
        case 'Method 3 (slowest)'
            methodid=2;
        otherwise
            return;
    end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    
    
    if isempty(sce.c_cell_type_tx) 
        if isempty(sce.c_cluster_id)
            errordlg('sce.c_cell_type_tx is empty.');
            return;        
        else
            answer = questdlg('sce.c_cell_type_tx is empty. Use sce.c_cluster_id?');
            if ~strcmp(answer,'Yes'), return; end
            cell_type_v=sce.c_cluster_id;                   
        end
    else
        answer = questdlg('Use sce.c_cell_type_tx?');
        switch answer
            case 'Yes'
                cell_type_v=sce.c_cell_type_tx;
            case 'No'
                if ~isempty(sce.c_cluster_id)
                    answer2 = questdlg('Use sce.c_cluster_id?');
                    if strcmp(answer2,'Yes')
                        cell_type_v=sce.c_cluster_id;                   
                    elseif strcmp(answer2,'No')
                        helpdlg('Action cancelled')
                        return;
                    else
                        return; 
                    end
                end
            otherwise
                return;
        end
    end
    
    
    [c,cL]=grp2idx(cell_type_v);
    if numel(cL)==1
        helpdlg('Only one cell type or cluster')
        return; 
    end    
    

    fw=gui.gui_waitbar;
    M=cell(numel(cL),2);
    [markerlist]=sc_pickmarkers(sce.X,sce.g,c,10,methodid);
    for k=1:numel(cL)        
        cLk=matlab.lang.makeValidName(cL{k});
        M{k,1}=cLk;
        M{k,2}=markerlist{k};
    end
    
% ==============    
    
    X=[]; szcl=[]; idcl=[];
    for k=1:length(cL)
        i=c==k;
        X=[X sce.X(:,i)];
        szcl=[szcl sum(i)];
        idcl=[idcl; c(i)];
    end
    X=sc_norm(X);
    X=log2(X+1);

% =========== 
Y=[]; idgn=[]; szgn=[]; Z=[];
% subi=1:10:size(X,2);
MX=[];
for k=1:numel(cL)    
    markerlist=M{k,2}(1:end);
    MX=[MX; markerlist];    
    [~,idx_g]=ismember(upper(markerlist),upper(sce.g));
    Y=[Y; X(idx_g,:)];
    idgn=[idgn; k*ones(length(markerlist),1)];
    szgn=[szgn length(markerlist)];
end


Y=zscore(Y,0,2);
qx=quantile(Y(:),0.90);
Y(Y>qx)=qx;
qx=quantile(Y(:),0.10);
Y(Y<qx)=qx;
    
    

Z=[];
for k=1:numel(cL)
    y=Y(idgn==k,:);
    
    z=[];
    for kk=1:numel(cL)
        z=[z mean(y(:,idcl==kk),2)];
    end
    %z1=grpstats(y.',idcl,@mean)';
    %assert(isequal(z,z1));
    Z=[Z; z];
end


    gui.gui_waitbar(fw);

f1=figure;
imagesc(Y);
szc=cumsum(szgn);
for k=1:max(idcl)-1
    xline(sum(idcl<k+1)+0.5,'r-');
    yline(szc(k)+0.5,'r-');
end
set(gca,'YTick',1:size(Y,1));
a=[]; b=[];
for k=1:max(idcl)
    a=[a sum(idcl<=k)];
    b=[b round(sum(idcl==k)./2)];
end
set(gca,'XTick',a-b);
set(gca,'XTickLabel',strrep(M(:,1),'_','\_'));
set(gca,'XTickLabelRotation',45);
set(gca,'YTick',1:length(MX));
set(gca,'YTickLabel',MX);
set(gca,'TickLength',[0 0])

tb1=uitoolbar(f1);
pt1 = uipushtool(tb1,'Separator','off');
pt1.Tooltip = 'Save marker gene map';
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','greencircleicon.gif'));
ptImage = ind2rgb(img,map);
pt1.CData = ptImage;
pt1.ClickedCallback = {@i_saveM,M};

pt1 = uipushtool(tb1,'Separator','off');
pt1.Tooltip = 'Summary map';
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object02.gif'));
ptImage = ind2rgb(img,map);
pt1.CData = ptImage;
pt1.ClickedCallback = @i_summarymap;

pt1 = uipushtool(tb1,'Separator','off');
pt1.Tooltip = 'Summary map, transposed';
[img,map] = imread(fullfile(matlabroot,...
            'toolbox','matlab','icons','HDF_object01.gif'));
ptImage = ind2rgb(img,map);
pt1.CData = ptImage;
pt1.ClickedCallback = @i_summarymapT;






    function i_summarymap(~,~)
        f2=figure;
        imagesc(Z);
        set(gca,'XTick',1:numel(cL));
        set(gca,'XTickLabel',strrep(M(:,1),'_','\_'));
        set(gca,'XTickLabelRotation',45);
        set(gca,'YTick',1:length(MX));
        set(gca,'YTickLabel',MX);
        set(gca,'TickLength',[0 0])
        f2.Position(1)=f2.Position(1)+200;
        f2.Position(2)=f2.Position(2)-200;
    end

    function i_summarymapT(~,~)
        f2=figure;
        imagesc(Z');
        set(gca,'YTick',1:numel(cL));
        set(gca,'YTickLabel',strrep(M(:,1),'_','\_'));
        set(gca,'XTickLabelRotation',45);
        set(gca,'XTick',1:length(MX));
        set(gca,'XTickLabel',MX);
        set(gca,'TickLength',[0 0])
        f2.Position(1)=f2.Position(1)+200;
        f2.Position(2)=f2.Position(2)-200;
    end

    function i_saveM(~,~,M)    
        labels = {'Save marker gene map M to variable named:'}; 
        vars = {'M'};
        values = {M};
        export2wsdlg(labels,vars,values);
    end    

end



%{
    FigureHandle=gcf;
    hAx = axes('Parent',FigureHandle);
    UitoolbarHandle = uitoolbar('Parent',FigureHandle ) ; 
    pt3 = uipushtool(UitoolbarHandle,'Separator','off');
    pt3.Tooltip = 'Select a gene to show expression';
    [img,map] = imread(fullfile(matlabroot,...
                'toolbox','matlab','icons','greencircleicon.gif'));
    ptImage = ind2rgb(img,map);
    pt3.CData = ptImage;
    pt3.ClickedCallback = {@i_saveM,M};

    function i_saveM(~,~,M)    
        labels = {'Save marker gene map M to variable named:'}; 
        vars = {'M'};
        values = {M};
        export2wsdlg(labels,vars,values);
    end
%}

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

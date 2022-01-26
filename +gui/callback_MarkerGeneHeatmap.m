function callback_MarkerGeneHeatmap(src,~)



    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    
%{    
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
%}

    [thisc,~]=gui.i_select1class(sce);
    if isempty(thisc)
        % errordlg('Undefined');
        return;
    end

    [c,cL]=grp2idx(thisc);
    if numel(cL)==1
        errordlg('Only one cell type or cluster.');
        return; 
    end    
    

    answer = questdlg('Generate marker gene heatmap',...
        'Select Method','Method 1 (DE 🐇)','Method 2 (scGeneFit 🐢)',...
        'Method 3 (LASSO 🐢🐢)','Method 1 (DE 🐇)');
    switch answer
        case 'Method 1 (DE 🐇)'
            methodid=1;
        case 'Method 2 (scGeneFit 🐢)'
            methodid=2;
        case 'Method 3 (LASSO 🐢🐢)'
            methodid=3;
        otherwise
            return;
    end    
    
    fw=gui.gui_waitbar;
    [markerlist]=sc_pickmarkers(sce.X,sce.g,c,10,methodid);
    M=cell(numel(cL),2);
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
    X=log(X+1);

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

% ======= customized heatmap - start
imagesc(Y);
szc=cumsum(szgn);
for k=1:max(idcl)-1
    xline(sum(idcl<k+1)+0.5,'r-');
    yline(szc(k)+0.5,'r-');
end
set(gca,'YTick',1:size(Y,1));
a=zeros(1,max(idcl)); b=zeros(1,max(idcl));
for k=1:max(idcl)
    a(k)=sum(idcl<=k);
    b(k)=round(sum(idcl==k)./2);
end
set(gca,'XTick',a-b);
% set(gca,'XTickLabel',strrep(M(:,1),'_','\_'));
set(gca,'XTickLabel',M(:,1));
set(gca,'XTickLabelRotation',45);
set(gca,'YTick',1:length(MX));
set(gca,'YTickLabel',MX);
set(gca,'TickLength',[0 0])
% ======= customized heatmap - end


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
        figure;
        h=heatmap(cL,MX,Z);
        h.Title = 'Marker Gene Heatmap';
        h.XLabel = 'Cell Type';
        h.YLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
    end

    function i_summarymapT(~,~)
         figure;
%         imagesc(Z');
%         set(gca,'YTick',1:numel(cL));
%         set(gca,'YTickLabel',strrep(M(:,1),'_','\_'));
%         set(gca,'XTickLabelRotation',45);
%         set(gca,'XTick',1:length(MX));
%         set(gca,'XTickLabel',MX);
%         set(gca,'TickLength',[0 0])
%         f2.Position(1)=f2.Position(1)+200;
%         f2.Position(2)=f2.Position(2)-200;
        h=heatmap(MX,cL,Z.');
        h.Title = 'Marker Gene Heatmap';
        h.YLabel = 'Cell Type';
        h.XLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
    end

    function i_saveM(~,~,M)
        if ~(ismcc || isdeployed)
            labels = {'Save marker gene map M to variable named:'}; 
            vars = {'M'};
            values = {M};
            export2wsdlg(labels,vars,values);
        else
            errordlg('This function is not available for standalone application.');
        end
    end    

end


function [cL]=e_sctype(sce)

[c,cL]=grp2idx(sce.c);
manuallyselect=false;
mfolder = fileparts(mfilename('fullpath'));
infile=fullfile(mfolder,'..','resources','ScTypeDB_full.xlsx');
if exist(infile,"file")
    T=readtable(infile);
end
utissuelist=unique(T.tissueType);
[indx1,tf1] = listdlg('PromptString',...
    {'Select Tissue Type(s):'},...
     'SelectionMode','multiple', ...
     'ListString',utissuelist,'ListSize',[220,300]);

if tf1~=1, return; end
selectedtissue=utissuelist(indx1);
Tm=T(ismember(T.tissueType,selectedtissue),:);
[Tw]=i_getTw(Tm);

% -------
wvalu=Tw.Var2;
wgene=string(Tw.Var1);
celltypev=string(Tm.cellName);
markergenev=string(Tm.geneSymbolmore1);

for k=1:max(c)
    ptsSelected=c==k;

    [Tct]=pkg.e_determinecelltype(sce, ptsSelected, wvalu, ...
        wgene, celltypev, markergenev);
    ctxt=Tct.C1_Cell_Type;

    if manuallyselect && length(ctxt)>1
        [indx, tf] = listdlg('PromptString', {'Select cell type'},...
            'SelectionMode', 'single', 'ListString', ctxt);
        if tf ~= 1, return; end
        ctxt = Tct.C1_Cell_Type{indx};
    else
        ctxt = Tct.C1_Cell_Type{1};
    end
    cL{k} = ctxt;
end

end


function [Tw]=i_getTw(Tm)
    s=upper(string(Tm.geneSymbolmore1));
    S=[];
    for k=1:length(s)
        a=strsplit(s(k),',');
        a=strtrim(a);    
        if strlength(a(end))==0 || isempty(a(end))
            a=a(1:end-1);
        end
        S=[S,a];
    end
    %%
    N=length(S);
    t=tabulate(S);
    f=cell2mat(t(:,3));
    if max(f)-min(f)<eps
        w=ones(N,1);
    else
        w=1+sqrt((max(f)-f)/(max(f)-min(f)));
    end
    genelist=string(t(:,1));
    Tw=table(genelist,w);
    Tw.Properties.VariableNames={'Var1','Var2'};    
end
%{
[Tm,Tw]=pkg.i_markerlist2weight(sce);
if isempty(Tm)||isempty(Tw)
    return;
end        
wvalu=Tw.Var2;
wgene=string(Tw.Var1);
celltypev=string(Tm.Var1);
markergenev=string(Tm.Var2);


[Tct]=i_determinecelltype(sce, ptsSelected, wvalu, ...
    wgene, celltypev, markergenev);


            ctxt=Tct.C1_Cell_Type;
            
            if manuallyselect && length(ctxt)>1
                [indx, tf] = listdlg('PromptString', {'Select cell type'},...
                    'SelectionMode', 'single', 'ListString', ctxt);
                if tf ~= 1, return; end
                ctxt = Tct.C1_Cell_Type{indx};
            else
                ctxt = Tct.C1_Cell_Type{1};
            end
            
            hold on;
            ctxtdisp = strrep(ctxt, '_', '\_');
            ctxtdisp = sprintf('%s_{%d}', ctxtdisp, i);
            cLdisp{i} = ctxtdisp;
            
            ctxt = sprintf('%s_{%d}', ctxt, i);
            cL{i} = ctxt;
            
            row = dataTipTextRow('', cLdisp(c));
            h.DataTipTemplate.DataTipRows = row;
            if size(sce.s, 2) >= 2
                siv = sce.s(ptsSelected, :);
                si = mean(siv, 1);
                idx = find(ptsSelected);
                [k] = dsearchn(siv, si);        % Nearest point search
                datatip(h, 'DataIndex', idx(k));
                % text(si(:,1),si(:,2),si(:,3),sprintf('%s',ctxt),...
                %     'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');
                %     elseif size(sce.s,2)==2
                %             si=mean(sce.s(ptsSelected,:));
                %             text(si(:,1),si(:,2),sprintf('%s',ctxt),...
                %                  'fontsize',10,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

%}
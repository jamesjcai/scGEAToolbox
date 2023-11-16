function callback_DPGene2Groups(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[i1, i2, cL1, cL2] = gui.i_select2grps(sce);
if length(i1) == 1 || length(i2) == 1, return; end


selitems={'TF Targets Expression',...    
    'MSigDB Signature',...
    'Predefined Signature'};
[indx1,tf1]=listdlg('PromptString',...
    'Select a metric for comparison.',...
    'SelectionMode','single','ListString',selitems, ...
    'ListSize',[200,300]);
if tf1~=1, return; end
[setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1);


[~,ix,iy]=intersect(upper(setgenes),upper(sce.g));
setmatrx=setmatrx(:,ix);    % s x g
X=sce.X(iy,:);              % g x c
Z=setmatrx*X;               % s x c


p_val=ones(size(Z,1),1);
avg_log2FC=nan(size(Z,1),1);
warning off
for k=1:size(Z,1)
    a=Z(k,i1);    
    b=Z(k,i2);
    p_val(k) = ranksum(a,b);
    if ~isnan(p_val(k))
        [ax]=nbinfit(a);
        [bx]=nbinfit(b);        
        avg_log2FC(k) = ax(1)./bx(1);    
    end
end
warning on
if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
end
T=table(setnames, avg_log2FC, p_val, p_val_adj);
T(isnan(T.p_val),:)=[];
T = sortrows(T, 'p_val_adj', 'ascend');
T=T(T.p_val_adj<0.01,:);

    outfile = sprintf('%s_vs_%s', ...
        matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));

    [~, filesaved] = gui.i_exporttable(T, true, 'T', outfile);    

    if ~isempty(filesaved)
        waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
    end
end


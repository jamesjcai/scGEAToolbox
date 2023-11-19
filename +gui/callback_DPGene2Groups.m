function callback_DPGene2Groups(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[i1, i2, cL1, cL2] = gui.i_select2grps(sce);
if length(i1) == 1 || length(i2) == 1, return; end



% selitems={'TF Targets Expression',...    
%     'MSigDB Signature',...
%     'Predefined Signature'};
% [indx1,tf1]=listdlg('PromptString',...
%     'Select a metric for comparison.',...
%     'SelectionMode','single','ListString',selitems, ...
%     'ListSize',[200,300]);
% if tf1~=1, return; end

[setmatrx, setnames, setgenes] = pkg.e_getgenesets; %(indx1);

fw = gui.gui_waitbar;

[~,ix,iy]=intersect(upper(setgenes),upper(sce.g));
setmatrx=setmatrx(:,ix);    % s x g
X=sce.X(iy,:);              % g x c
Z=setmatrx*X;               % s x c
gsetsize=sum(setmatrx,2);   % gene number per set

p_val=ones(size(Z,1),1);
avg_log2FC=nan(size(Z,1),1);
v1=nan(size(Z,1),1);
v2=nan(size(Z,1),1);
n1=nan(size(Z,1),1);
n2=nan(size(Z,1),1);
m1=nan(size(Z,1),1);
m2=nan(size(Z,1),1);

warning off
for k=1:size(Z,1)
    a=Z(k,i1);    
    b=Z(k,i2);
    p_val(k) = ranksum(a,b);
    if ~isnan(p_val(k)) && p_val(k)<1e-3        
        %[ax]=nbinfit(a);
        %[bx]=nbinfit(b);

        [ax]=mean(a);
        [bx]=mean(b);
        avg_log2FC(k) = log2(ax(1)./bx(1));
        v1(k)=ax(1);
        v2(k)=bx(1);
        n1(k)=numel(a);
        n2(k)=numel(b);
        m1(k)=sum(a>0);
        m2(k)=sum(b>0);
    end
end
warning on
if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
end

T=table(setnames, gsetsize, v1, v2, avg_log2FC, m1, n1, m2, n2, p_val, p_val_adj);
T(isnan(T.p_val)|isnan(T.avg_log2FC)|T.avg_log2FC<1,:)=[];
T = sortrows(T, 'p_val_adj', 'ascend');
T=T(T.p_val_adj<0.01 & T.gsetsize>5,:);

    gui.gui_waitbar(fw);

    if size(T,1)==0
        waitfor(helpdlg('No significant results.',''));
    else
        outfile = sprintf('%s_vs_%s', ...
            matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
        [~, filesaved] = gui.i_exporttable(T, true, 'T', outfile);    
        if ~isempty(filesaved)
            waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
        end
    end
end


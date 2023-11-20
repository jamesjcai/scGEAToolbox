function callback_DPGene2Groups(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[i1, i2, cL1, cL2] = gui.i_select2grps(sce);
if length(i1) == 1 || length(i2) == 1, return; end

c=zeros(size(i1));
c(i1)=1; c(i2)=2;
cL=[cL1;cL2];
if ~all(c>0)
    sce=sce.selectcells(c>0);
    c=c(c>0);
    i1=c==1;
    i2=c==2;
end


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

[gc,ix,iy]=intersect(upper(setgenes),upper(sce.g)); %,'stable');
setgenes=setgenes(ix);
setmatrx=setmatrx(:,ix);    % s x g

X=sce.X(iy,:);              % g x c
sce.X=X;
sce.g=sce.g(iy);

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

%warning off
for k=1:size(Z,1)
    if any(setmatrx(k,:))
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
end
%warning on
if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
end

T=table(setnames, gsetsize, v1, v2, avg_log2FC, m1, n1, m2, n2, p_val, p_val_adj);
T(isnan(T.p_val)|isnan(T.avg_log2FC)|abs(T.avg_log2FC)<1,:)=[];
T = sortrows(T, 'p_val_adj', 'ascend');
T=T(T.p_val_adj<0.01 & T.gsetsize>=5,:);

    gui.gui_waitbar(fw);

    if size(T,1)==0
        waitfor(helpdlg('No significant results.',''));
        return;
    else
        outfile = sprintf('%s_vs_%s', ...
            matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
        %filesaved = fullfile(outdir, outfile);
        %writetable(T, filesaved, 'FileType', 'spreadsheet');
        [~, filesaved] = gui.i_exporttable(T, true, 'T', outfile);    
        if ~isempty(filesaved)
           waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
        end
    end



answer=questdlg('Plot results? You will be asked to select a folder to save figure files. Continue?','');
if ~strcmp(answer,'Yes'), return; end
outdir = uigetdir;
if ~isfolder(outdir), return; end
    


%                c=zeros(size(i1));
%                c(i1)=1; c(i2)=2; 
%                sce=sce.selectcells(c>0);

                
                Xt=log(1+sc_norm(sce.X));

 fw = gui.gui_waitbar_adv;
 for k=1:size(T,1)
     
     idx=T.setnames(k)==setnames;
     posg=string(setgenes(setmatrx(idx,:)));
    
% % [y] = gui.e_cellscore(sce, posg);
% % [i1, i2, cL1, cL2] = gui.i_select2grps(sce);
                % gui.i_violinplot(y, thisc, ttxt, true, [], posg);
                % xlabel('Cell group');
                % ylabel('Cellular score');
                
                %figure;
        outfile1 = sprintf('dotplot_%s.png', ...
            matlab.lang.makeValidName(T.setnames(k)));
        outfile2 = sprintf('violplt_%s.png', ...
            matlab.lang.makeValidName(T.setnames(k)));
        
        filesaved1 = fullfile(outdir, outfile1);
        filesaved2 = fullfile(outdir, outfile2);

        assignin("base","Xt",Xt);
        assignin("base","g",sce.g);
        assignin("base","c",c);
        assignin("base","cL",cL);
        assignin("base","posg",posg);
        gui.gui_waitbar_adv(fw,(k-1)./size(T,1));

                %try
                    f1=gui.i_dotplot(Xt, upper(sce.g), c, cL, upper(posg), true, T.setnames(k));
                    saveas(f1, filesaved1);
                %catch ME
                %    warning(ME.message);
                %end
                
                %try
                    [y] = gui.e_cellscore(sce, posg, 2, false);  % 'AddModuleScore/Seurat'
                    ttxt=T.setnames(k);
                    f2=gui.i_violinplot(y, c, ttxt, true, [], posg);
                    saveas(f2, filesaved2);
                %catch ME
                %    warning(ME.message);
                %end
                    

 end
 gui.gui_waitbar_adv(fw);

end


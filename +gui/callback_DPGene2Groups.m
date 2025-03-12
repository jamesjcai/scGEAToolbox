function callback_DPGene2Groups(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('DP Analysis', FigureHandle), return; end

    extprogname = 'scgeatool_DPAnalysis';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end


[i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2), return; end

% --------
a=sprintf('%s vs. %s',cL1{1}, cL2{1});
b=sprintf('%s vs. %s',cL2{1}, cL1{1});
answer = gui.myQuestdlg(FigureHandle, 'Which vs. which?','',{a,b},a);
switch answer
    case a
    case b
        i3=i1; i1=i2; i2=i3;
        cL3=cL1; cL1=cL2; cL2=cL3;
    otherwise
        return;
end
% ----------


c=zeros(size(i1));
c(i1)=1; c(i2)=2;
cL=[cL1;cL2];
if ~all(c>0)
    sce = sce.selectcells(c>0);
    c = c(c>0);
    i1 = c==1;
    i2 = c==2;
end

[indx1,speciesid] = gui.i_selgenecollection;
if isempty(indx1), return; end
[setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1, speciesid); %(indx1);
if isempty(setmatrx) || isempty(setnames) || isempty(setgenes) 
    return;
end

assignin("base", "setmatrx", setmatrx);
assignin("base", "setnames", setnames);
assignin("base", "setgenes", setgenes);

fw = gui.myWaitbar(FigureHandle);

[~,ix,iy]=intersect(upper(setgenes), ...
                    upper(sce.g)); %,'stable');
setgenes=setgenes(ix);
setmatrx=setmatrx(:,ix);      % s x g

sceX = log1p(sc_norm(sce.X));
X = sceX(iy,:);              % g x c


sce.X = sce.X(iy,:);
sce.g = sce.g(iy);

Z = setmatrx*X;               % s x c
gsetsize = sum(setmatrx,2);   % gene number per set

p_val = ones(size(Z,1),1);
avg_log2FC = nan(size(Z,1),1);
v1 = nan(size(Z,1),1);
v2 = nan(size(Z,1),1);
n1 = nan(size(Z,1),1);
n2 = nan(size(Z,1),1);
m1 = nan(size(Z,1),1);
m2 = nan(size(Z,1),1);

%warning off
for k = 1:size(Z,1)
    if any(setmatrx(k,:))
        a = Z(k,i1);
        b = Z(k,i2);
        p_val(k) = ranksum(a,b);
        if ~isnan(p_val(k)) && p_val(k)<1e-3
            %[ax]=nbinfit(a);
            %[bx]=nbinfit(b);
            [ax] = mean(a);
            [bx] = mean(b);
            avg_log2FC(k) = log2(ax(1)./bx(1));
            v1(k) = ax(1);
            v2(k) = bx(1);
            n1(k) = numel(a);
            n2(k) = numel(b);
            m1(k) = sum(a>0);
            m2(k) = sum(b>0);
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

    gui.myWaitbar(FigureHandle, fw);


    if height(T)==0
        gui.myHelpdlg(FigureHandle, 'No significant results.');
        return;
    else
        outfile = sprintf('%s_vs_%s_DP_results', ...
            matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
        %filesaved = fullfile(outdir, outfile);
        %writetable(T, filesaved, 'FileType', 'spreadsheet');
        [~, filesaved] = gui.i_exporttable(T, true, 'Tdpgenelist', outfile);
            % "Tcellattrib","CellAttribTable"
            % "Tviolindata","ViolinPlotTable"
            % "Tcrosstabul","CrosstabulTable"
            % "Tcellsignmt","CellSignatTable"
            % "Tdpgenesres","DPGenesResTable"
        if ~isempty(filesaved)
           gui.myHelpdlg(FigureHandle, sprintf('Result has been saved in %s',filesaved));
        end
    end

answer=gui.myQuestdlg(FigureHandle, 'Select gene sets and plot results?','');
if ~strcmp(answer,'Yes')
    return;
else
    [idxneedplot] = gui.i_selmultidlg(T.setnames, [], FigureHandle);
end
if isempty(idxneedplot), return; end

outdir = tempdir;

Xt=log1p(sc_norm(sce.X));
images = {};

 fw = gui.gui_waitbar_adv;
 success=false;

 for kk=1:numel(idxneedplot)
     k=idxneedplot(kk);
     if ~ismember(k,idxneedplot), continue; end
     idx=T.setnames(k)==setnames;
     posg=string(setgenes(setmatrx(idx,:)));
    
        outfile1 = sprintf('dotplot_%s.png', ...
            matlab.lang.makeValidName(T.setnames(k)));
        outfile2 = sprintf('violplt_%s.png', ...
            matlab.lang.makeValidName(T.setnames(k)));
        
        filesaved1 = fullfile(outdir, outfile1);
        filesaved2 = fullfile(outdir, outfile2);

        % assignin("base","Xt",Xt);
        % assignin("base","g",sce.g);
        % assignin("base","c",c);
        % assignin("base","cL",cL);
        % assignin("base","posg",posg);
        
        % gui.gui_waitbar_adv(fw,(k-1)./height(T));
        gui.gui_waitbar_adv(fw,(kk-1)./numel(idxneedplot));
        
    suc1=false;
    try
        f1 = gui.i_dotplot(Xt, upper(sce.g), c, cL, upper(posg), true, T.setnames(k));
        saveas(f1, filesaved1);                    
        images = [images, {filesaved1}];
        suc1=true;
    catch ME
        warning(ME.message);
    end
    
    suc2=false;
    try
        [y] = gui.e_cellscore(sce, posg, 2, false);  % 'AddModuleScore/Seurat'
        ttxt=T.setnames(k);
        % assignin("base","y",y);
        % assignin("base","c",c);
        % assignin("base","cL",cL);
        % assignin("base","posg",posg);

        f2 = gui.i_violinplot(y, cL(c), ttxt, true, [], posg);
        saveas(f2, filesaved2);
        images = [images, {filesaved2}];
        suc2 = true;
    catch ME
        %success=false;
        warning(ME.message);
    end

    success = suc1 & suc2;
 end
 gui.gui_waitbar_adv(fw);
 if ~success
     gui.myHelpdlg(FigureHandle, 'All figure files are not saved.');
     winopen(outdir);
 end

    
    % answer = gui.myQuestdlg(FigureHandle, 'Output to PowerPoint?','','Yes','No','Yes');
    % switch answer
    %     case 'Yes'
    %         needpptx = true;
    %     case 'No'
    %         needpptx = false;
    %     otherwise
    %         needpptx = false;
    % end

    % if needpptx
        gui.i_save2pptx(images);
    % else
        % if success    
        %     answer = gui.myQuestdlg(FigureHandle, sprintf('Figure files have been saved in %s. Open the folder to view files?', outdir),'');
        %     if strcmp(answer, 'Yes'), winopen(outdir); end
        % end
    % end
end

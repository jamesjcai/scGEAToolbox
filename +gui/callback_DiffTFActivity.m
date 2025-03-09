function callback_DiffTFActivity(src, ~)

[~, sce, isui] = gui.gui_getfigsce(src);

[thisc, clabel] = gui.i_select1class(sce, false);
if isempty(thisc), return; end
species = gui.i_selectspecies(2);
if isempty(species), return; end

answer = questdlg('Select algorithm:', '', ...
    'UCell [PMID:34285779]', 'NMF [PMID:33135076]🐢', ...
    'UCell [PMID:34285779]');
if isempty(answer), return; end

switch answer
    case 'NMF [PMID:33135076]🐢'
        methodid = 3;
    case 'UCell [PMID:34285779]'
        methodid = 1;
    otherwise
        return;
end


if methodid ~= 3
    fw = gui.gui_waitbar;
end

[cs, tflist, ~, numtargetgenes] = sc_tfactivity(sce.X, sce.g, [], species, methodid);
T = pkg.e_grptest(cs, thisc);
T = [table(tflist), T, table(numtargetgenes)];

if methodid ~= 3
    gui.gui_waitbar(fw);
end

upperg = upper(sce.g);

[tfexpressed] = ismember(upper(tflist), upperg);
T = [T, table(tfexpressed)];
T2 = T(tfexpressed, :);
cs2 = cs(tfexpressed, :);

outfile = sprintf('DiffTFActivity_%s', matlab.lang.makeValidName(clabel));

answer = questdlg(sprintf('%d out of %d TFs expressed in cells. Keep only %d expressed TFs?', ...
    height(T2), height(T), height(T2)));
switch answer
    case 'Yes'
        T = T2;
        cs = cs2;
    case 'No'
    otherwise
        return;
end


tflist = string(T.tflist);
regdirection = zeros(length(tflist), 1);
for k = 1:length(tflist)
    [tfexpressed, idx] = ismember(upper(tflist(k)), upperg);
    if tfexpressed
        thiscs = cs(k, :);
        thisex = sce.X(idx, :);
        regdirection(k) = corr(thisex(:), thiscs(:), 'type', 'Spearman');
    end
end
T = [T, table(regdirection)];


try
    if length(unique(thisc)) > 2
        [T, idx] = sortrows(T, "p_anova", "ascend");
        cs = cs(idx, :);
        [T, idx] = sortrows(T, "p_kruskalwallis", "ascend");
        cs = cs(idx, :);
    else
        [T, idx] = sortrows(T, "p_ttest", "ascend");
        cs = cs(idx, :);
        [T, idx] = sortrows(T, "p_wilcoxon", "ascend");
        cs = cs(idx, :);
    end
catch

end
gui.i_exporttable(T, true, 'T', outfile);

answer = questdlg('Violin plot for top TFs with most variable activity between groups?');

    switch answer
        case 'Yes'
            [numfig] = gui.i_inputnumg(length(T.tflist));
            if isempty(numfig) || isnan(numfig), return; end
            if isnumeric(numfig) && numfig > 0 && numfig <= length(T.tflist)

                %[~,cL]=grp2idx(thisc);
                [~, cL, noanswer] = gui.i_reordergroups(thisc);
                if noanswer, return; end
                %                 cL=strrep(cL,'_','\_');
                %                 thisc=strrep(thisc,'_','\_');

                F = cell(numfig, 1);
                for k = 1:numfig
                    ttxt = T.tflist(k);
                    f = gui.i_violinplot(cs(k, :), thisc, ttxt, true, cL);
                    %                     f = figure('visible','off');
                    %                     pkg.i_violinplot(cs(k,:),thisc,true,cL);
                    %                     title(T.tflist(k));
                    %                     ylabel('TF activity')
                    %                     P = get(f,'Position');
                    %                     set(f,'Position',[P(1)-20*k P(2)-20*k P(3) P(4)]);
                    %                     set(f,'visible','on');
                    %                     %f.Position(3) = f.Position(3) * 2.2;
                    %                     drawnow;
                    F{k} = f;
                end
                gui.i_export2pptx(F, string(T.tflist(1:numfig)));
            end
    end
end

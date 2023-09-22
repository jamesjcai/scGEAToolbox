function [f] = i_violinplot(y, thisc, ttxt, colorit, cLorder, posg)

if nargin < 6, posg = []; end % used in callback_CompareGeneBtwCls
if nargin < 5 || isempty(cLorder)
    [~, cLorder] = grp2idx(thisc);
end
if nargin < 4
    colorit = true;
end
f = figure('visible', 'off');
tb = uitoolbar(f);
pkg.i_addbutton2fig(tb, 'off', {@i_savedata, y, thisc}, ...
    'export.gif', 'Export data...');
pkg.i_addbutton2fig(tb, 'off', {@i_testdata, y, thisc}, ...
    'icon-fa-stack-exchange-10.gif', 'ANOVA/T-test...');
pkg.i_addbutton2fig(tb, 'off', @i_addsamplesize, ...
    "icon-mat-blur-linear-10.gif", 'Add Sample Size');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, ...
    "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', @i_invertcolor, ...
    "plotpicker-pie.gif", 'Switch BW/Color');
pkg.i_addbutton2fig(tb, 'off', @i_reordersamples, ...
    "plotpicker-errorbar.gif", 'Reorder Samples');
pkg.i_addbutton2fig(tb, 'off', @i_sortbymean, ...
    "plotpicker-cra.gif", 'Sort Samples by Median');
pkg.i_addbutton2fig(tb, 'off', @i_renametitle, ...
    "icon-mat-touch-app-10.gif", 'Change Plot Title');
pkg.i_addbutton2fig(tb, 'on', @i_viewgenenames, ...
    'HDF_point.gif', 'Show Gene Names');


isdescend = false;
OldTitle = [];
OldXTickLabel = [];
cLorder = strrep(cLorder, '_', '\_');
thisc = strrep(string(thisc), '_', '\_');
pkg.i_violinplot(y, thisc, colorit, cLorder);
title(strrep(ttxt, '_', '\_'));
%ylabel(selitems{indx1});
movegui(f, 'center');
set(f, 'visible', 'on');

%catch ME
%    errordlg(ME.message);
%end

    function i_renametitle(~, ~)
        helpdlg('Double-click on the title to make change.', '');
end

        function i_invertcolor(~, ~)
            colorit = ~colorit;
            cla;
            pkg.i_violinplot(y, thisc, colorit, cLorder);
    end

            function i_addsamplesize(~, ~)
                b = gca;
                if isempty(OldXTickLabel)
                    a = zeros(length(cLorder), 1);

                    OldXTickLabel = b.XTickLabel;
                    for k = 1:length(cLorder)
                        a(k) = sum(thisc == cLorder(k));
                        b.XTickLabel{k} = sprintf('%s\\newline(n=%d)', ...
                            b.XTickLabel{k}, a(k));
                    end
                else
                    b.XTickLabel = OldXTickLabel;
                    OldXTickLabel = [];
                end
        end

                function i_sortbymean(~, ~)
                    [cx, cLx] = grp2idx(thisc);
                    a = zeros(max(cx), 1);
                    for k = 1:max(cx)
                        a(k) = median(y(cx == k));
                    end
                    if isdescend
                        [~, idx] = sort(a, 'ascend');
                        isdescend = false;
                    else
                        [~, idx] = sort(a, 'descend');
                        isdescend = true;
                    end
                    cLx_sorted = cLx(idx);

                    %[~,cL,noanswer]=gui.i_reordergroups(thisc,cLx_sorted);
                    %if noanswer, return; end
                    cla
                    cLorder = cLx_sorted;
                    pkg.i_violinplot(y, thisc, colorit, cLorder);
            end


                    function i_reordersamples(~, ~)
                        [~, cLorder, noanswer] = gui.i_reordergroups(thisc);
                        if noanswer, return; end
                        cla
                        pkg.i_violinplot(y, thisc, colorit, cLorder);
                end


                        function i_testdata(~, ~, y, grp)
                            a = gca;
                            if isempty(OldTitle)
                                OldTitle = a.Title.String;
                                if size(y, 2) ~= length(grp)
                                    y = y.';
                                end
                                tbl = pkg.e_grptest(y, grp);
                                %h1=gca;
                                %titre=string(h1.Title.String);

                                %     a=sprintf('%s\n%s=%.2e; %s=%.2e', ...
                                %         strrep(string(ttxt),'_','\_'), ...
                                %         strrep(tbl.Properties.VariableNames{1},'_','\_'), ...
                                %         tbl.(tbl.Properties.VariableNames{1}), ...
                                %         strrep(tbl.Properties.VariableNames{2},'_','\_'), ...
                                %         tbl.(tbl.Properties.VariableNames{2}));
                                %     title(a);

                                b = sprintf('%s=%.2e; %s=%.2e', ...
                                    strrep(tbl.Properties.VariableNames{1}, '_', '\_'), ...
                                    tbl.(tbl.Properties.VariableNames{1}), ...
                                    strrep(tbl.Properties.VariableNames{2}, '_', '\_'), ...
                                    tbl.(tbl.Properties.VariableNames{2}));
                                if iscell(OldTitle)
                                    newtitle = OldTitle;
                                else
                                    newtitle = {OldTitle};
                                end
                                newtitle{2} = b;
                                a.Title.String = newtitle;
                            else
                                a.Title.String = OldTitle;
                                OldTitle = [];
                                %gui.i_exporttable(tbl,true);
                            end
                    end


                            function i_viewgenenames(~, ~)
                                if isempty(posg)
                                    helpdlg('The gene set is empty. This score may not be associated with any gene set.');
                                else
                                    %idx=matches(sce.g,posg,'IgnoreCase',true);
                                    %gg=sce.g(idx);
                                    inputdlg(ttxt, ...
                                        '', [10, 50], ...
                                        {char(posg)});
                                end
                        end

                        end

                            function i_savedata(~, ~, a, b)
                            T = table(a(:), b(:));
                            T.Properties.VariableNames = {'ScoreLevel', 'GroupID'};
                            %T=sortrows(T,'ScoreLevel','descend');
                            %T=sortrows(T,'GroupID');
                            gui.i_exporttable(T, true);
                        end

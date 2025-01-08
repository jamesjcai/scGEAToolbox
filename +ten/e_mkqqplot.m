function e_mkqqplot(T)

pd = makedist('Gamma', 'a', 0.5, 'b', 2);

hx=gui.myFigure;

hFig = hx.FigureHandle;
a = gca(hFig);
qqplot(a, T.FC, pd);

[~, idx] = sort(T.FC);
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1x, T.genelist(idx)};


% tb = findall(hFig, 'Tag', 'FigureToolBar'); 
%tb = uitoolbar('Parent', hFig);
%uipushtool(tb, 'Separator', 'off');

%tb = uitoolbar('Parent', hFig);
%hx.addCustomButton('off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
%hx.addCustomButton('on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
%hx.addCustomButton('off', {@gui.i_savemainfig, 2}, "svg.gif", 'Save Figure as Graphic File...');
%hx.addCustomButton('off', {@gui.i_savemainfig, 1}, "svg.gif", 'Save Figure as SVG File...');
hx.addCustomButton('on', @i_savetable, 'export.gif', 'Export data...');
%hx.addCustomButton('on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');
hx.show();


% h1=h(1);
% h1.DataTipTemplate.DataTipRows = T.genelist(idx);
% for k=1:5
%     datatip(h1, 'DataIndex', idx(k));
% end

    function i_savetable(~, ~)
        answer = questdlg('Export & save data to:', '', ...
            'Workspace', 'TXT/CSV file', 'Excel file', 'Workspace');
        if ~isempty(answer)
            switch answer
                case 'Workspace'
                    labels = {'Save T to variable named:'};
                    vars = {'T'};
                    values = {T};
                    [~, ~] = export2wsdlg(labels, vars, values, ...
                        'Save Data to Workspace');
                case 'TXT/CSV file'
                    [file, path] = uiputfile({'*.csv'; '*.*'}, 'Save as');
                    if isequal(file, 0) || isequal(path, 0)
                        return;
                    else
                        fw = gui.gui_waitbar;
                        filename = fullfile(path, file);
                        writetable(T, filename, 'FileType', 'text');
                        gui.gui_waitbar(fw);
                    end
                case 'Excel file'

                    [file, path] = uiputfile({'*.xlsx'; '*.*'}, 'Save as');
                    if isequal(file, 0) || isequal(path, 0)
                        return;
                    else
                        fw = gui.gui_waitbar;
                        filename = fullfile(path, file);
                        writetable(T, filename, 'FileType', 'spreadsheet');
                        gui.gui_waitbar(fw);
                    end
            end
        end
    end



end

function txt = i_myupdatefcn1x(~, event_obj, g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end

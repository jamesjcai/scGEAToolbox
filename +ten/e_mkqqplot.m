function e_mkqqplot(T)

pd = makedist('Gamma', 'a', 0.5, 'b', 2);

hx=gui.myFigure;

hFig = hx.FigureHandle;
a = gca(hFig);
qqplot(a, T.FC, pd);

[~, idx] = sort(T.FC);
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1x, T.genelist(idx)};
hx.addCustomButton('off', @i_savetable, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.show();


% h1=h(1);
% h1.DataTipTemplate.DataTipRows = T.genelist(idx);
% for k=1:5
%     datatip(h1, 'DataIndex', idx(k));
% end

    function i_savetable(~, ~)
        answer = gui.myQuestdlg(hFig, 'Export & save data to:', '', ...
            {'Workspace', 'TXT/CSV file', 'Excel file'}, 'Workspace');
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

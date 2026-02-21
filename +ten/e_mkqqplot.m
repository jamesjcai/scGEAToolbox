function e_mkqqplot(T)


% T = T(2:end,:);
pd = makedist('Gamma', 'a', 0.5, 'b', 2);

hx=gui.myFigure;
hFig = hx.FigHandle;
if gui.i_isuifig(hFig)
    a = hx.AxHandle;
else
    a = gca;
end


h = qqplot(a, T.FC, pd);
[~, idx] = sort(T.FC);
sorted_g = T.genelist(idx);

dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn1x, sorted_g};
hx.addCustomButton('off', @in_callback_savetable, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.show();

assignin("base","h",h)


 
h1=h(1);
%h.DataTipTemplate.DataTipRows = T.genelist(idx);
for k=0:5
    %T.genelist{end-k}
    %h1.XData(end-k)
    %h1.YData(end-k)
     % datatip(h, 'DataIndex', idx(k));
     text(h1.XData(end-k), h1.YData(end-k), sorted_g{end-k}, 'Rotation', 45);
     % annotation("textarrow", h1.XData(end-k), h1.YData(end-k),'String', "aa");
end



    function in_callback_savetable(~, ~)
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
if iscell(g(idx))
    txt = g(idx);
else
    txt = {g(idx)};
end
end

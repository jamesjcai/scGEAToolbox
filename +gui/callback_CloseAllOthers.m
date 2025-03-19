function callback_CloseAllOthers(src, ~)
    % Callback function that closes all other MATLAB figures except the current one.
    
    % Get the handle of the currently active figure
    % FigureHandle = gcf;
    [FigureHandle, ~, isui] = gui.gui_getfigsce(src);
    
    % Find all open figure handles in the workspace
    allFigures = findobj(0, 'type', 'figure');
    
    % Check if the current figure is in the list of all figures and there's more than one figure
    [isCurrentInList, currentIndex] = ismember(FigureHandle, allFigures);
    
    if isCurrentInList && length(allFigures) > 1
        % Ask the user to confirm before closing other figures
        confirmation = gui.myQuestdlg(FigureHandle, 'Close all other figures?', 'Confirmation');
        if isempty(confirmation), return; end
        % If the user does not confirm, exit the function
        if ~strcmp(confirmation, 'Yes')
            return;
        end
        
        % Loop through all figure handles and close the ones that are not the current one
        for index = 1:length(allFigures)
            if index ~= currentIndex
                try
                    % Attempt to close the figure
                    close(allFigures(index));
                    
                catch closeError
                    % Handle any exceptions that occur during closing a figure
                    disp(['Failed to close figure ', num2str(allFigures(index)), ': ', closeError.message]);
                end
            end
        end
        if isui
            uialert(FigureHandle, 'All other figures have been closed.','','Icon','info');
        else
            disp('All other figures have been closed.');
        end
    else
         if isui
            uialert(FigureHandle, 'All other figures have been closed.','','Icon','info');
        else
            disp('Either there is only one figure or the current figure is part of another application group. Nothing has been done.');
         end
         
    end
end

%{
function callback_CloseAllOthers(~, ~)
a = gcf;
all_figs = findobj(0, 'type', 'figure');
[y, idx] = ismember(a, all_figs);
if y && length(all_figs) > 1
    answer = gui.myQuestdlg(FigureHandle, 'Close all other figures?');
    if ~strcmp(answer, 'Yes'), return; end
    for k = 1:length(all_figs)
        if k ~= idx
            try
                close(all_figs(k));
            catch
            end
        end
    end
end
end
%}
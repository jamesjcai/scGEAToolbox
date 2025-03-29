function callback_CloseAllOthers(src, ~)
    % Callback function that closes all other MATLAB figures except the current one.
    
    % Get the handle of the currently active figure
    [FigureHandle] = gui.gui_getfigsce(src);
    
    allFigures = findall(0, 'Type', 'Figure');

    % Check if the current figure is in the list of all figures and there's more than one figure
    [isCurrentInList, currentIndex] = ismember(FigureHandle, allFigures);
    
    if isCurrentInList && length(allFigures) > 1
        % Ask the user to confirm before closing other figures
        confirmation = gui.myQuestdlg(FigureHandle, ...
            'Close all other figures?', 'Confirmation');
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
        gui.myHelpdlg(FigureHandle, 'All other figures have been closed.', '', true);
    end
end


function i_renametitle(~, ~)
    ax = gca; % Get the current axes
    titleObj = ax.Title; % Access the Title property of the axes
    
    if isempty(titleObj.String)
        title(ax, 'Title');
        % %fprintf('The current plot does not have a title.\n');
        % prompt = {'Enter a title for the plot:'};
        % dlgTitle = 'Set Plot Title';
        % dims = [1 50]; % Single-line input dialog
        % defaultAnswer = {''};
        % userInput = inputdlg(prompt, dlgTitle, dims, defaultAnswer);
        % 
        % % Check if user entered a title
        % if ~isempty(userInput) && ~isempty(userInput{1})
        %     newTitle = userInput{1}; % Extract the entered title
        %     title(newTitle); % Set the new title
        %     %fprintf('The title has been set to: "%s"\n', newTitle);
        % else
        %     %fprintf('No title was set since the input was empty.\n');
        % end
    end
    helpdlg('Double-click on the title to make change.', '');
end

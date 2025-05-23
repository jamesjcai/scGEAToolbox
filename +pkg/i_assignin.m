function saveName = i_assignin(data, baseName)
    % SAVETOWORKSPACEWITHAUTONAME Saves a variable to the base workspace
    % with automatic renaming if the name already exists
    %
    % Inputs:
    %   data - The data to save to the workspace
    %   baseName - The desired variable name (string)
    %
    % Output:
    %   saveName - The actual name used for saving (may be renamed)
    
    % Check if base name exists in workspace
    if evalin('base', sprintf('exist(''%s'', ''var'')', baseName))
        % Name exists, find a new name
        counter = 1;
        while true
            newName = sprintf('%s_%d', baseName, counter);
            if ~evalin('base', sprintf('exist(''%s'', ''var'')', newName))
                % Found an available name
                saveName = newName;
                break;
            end
            counter = counter + 1;
        end
    else
        % Base name is available
        saveName = baseName;
    end
    
    % Save the variable to workspace with the determined name
    assignin('base', saveName, data);
    
    % Output message (optional)
    fprintf('Variable saved as: %s\n', saveName);
end
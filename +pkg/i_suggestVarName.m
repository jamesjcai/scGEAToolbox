function newVarName = i_suggestVarName(baseName, workspace)
    % SUGGESTVARNAME Suggests a new variable name by appending numbers
    % if the base name already exists in the workspace
    %
    % Usage:
    %   newVarName = suggestVarName('myVar')           % Uses caller workspace
    %   newVarName = suggestVarName('myVar', 'base')   % Uses base workspace
    %   newVarName = suggestVarName('myVar', 'caller') % Uses caller workspace
    %
    % Input:
    %   baseName  - string or char array with the desired base variable name
    %   workspace - (optional) 'base', 'caller', or 'auto' (default: 'auto')
    %               'auto' uses base workspace if called from command line,
    %               caller workspace if called from function
    %
    % Output:
    %   newVarName - string with suggested variable name

    % Handle input arguments
    if nargin < 2
        workspace = 'auto';
    end

    % Convert to string if input is char array
    if ischar(baseName)
        baseName = string(baseName);
    end

    % Determine which workspace to use
    if strcmp(workspace, 'auto')
        % Check if we're in the command window or a function
        try
            % Try to get caller workspace variables
            evalin('caller', 'who');
            workspace = 'caller';
        catch
            % If no caller workspace, use base workspace
            workspace = 'base';
        end
    end

    % Get all existing variable names from the target workspace
    try
        switch lower(workspace)
            case 'base'
                existingVars = evalin('base', 'who');
            case 'caller'
                existingVars = evalin('caller', 'who');
            otherwise
                error('Invalid workspace. Use: base, caller, or auto');
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:evalin:invalidWorkspace')
            % Fallback to base workspace if caller doesn't exist
            existingVars = evalin('base', 'who');
        else
            rethrow(ME);
        end
    end

    % Convert existing variable names to strings for comparison
    existingVars = string(existingVars);

    % Check if base name exists
    if ~any(strcmp(existingVars, baseName))
        % Base name doesn't exist, return it as is
        newVarName = baseName;
        return;
    end

    % Base name exists, find next available name by checking all existing variables
    counter = 1;
    while true
        candidateName = baseName + string(counter);

        % Check if this candidate exists in the variable list
        if ~any(strcmp(existingVars, candidateName))
            newVarName = candidateName;
            return;
        end

        counter = counter + 1;

        % Safety check to prevent infinite loops
        if counter > 10000
            error('Too many existing variables with similar names (>10000)');
        end
    end
end

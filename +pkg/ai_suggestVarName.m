function newVarName = ai_suggestVarName(baseName, workspace)
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

% Example usage and demonstration
function demonstrateVarSuggestion()
    % Create some variables in workspace
    myVar = 1;
    myVar1 = 2;
    myVar2 = 3;
    data = [1, 2, 3];
    data3 = [4, 5, 6];  % Note: data1 and data2 don't exist, but data3 does

    % Get current workspace variables
    currentVars = who;
    fprintf('Current workspace variables: %s\n', strjoin(currentVars, ', '));

    % Test the function with different scenarios
    fprintf('\n=== Testing variable name suggestions ===\n');

    % Test 1: Variable exists, should suggest next number
    suggested1 = suggestVarName('myVar');
    fprintf('Suggested name for "myVar": %s\n', suggested1);

    % Test 2: Variable doesn't exist, should return as is
    suggested2 = suggestVarName('newVar');
    fprintf('Suggested name for "newVar": %s\n', suggested2);

    % Test 3: Gap in numbering (data1, data2 don't exist but data3 does)
    suggested3 = suggestVarName('data');
    fprintf('Suggested name for "data": %s\n', suggested3);

    % Test 4: Test with specific workspace
    suggested4 = suggestVarName('testVar', 'base');
    fprintf('Suggested name for "testVar" in base workspace: %s\n', suggested4);

    % Test 5: Create variable using suggested name
    suggestedName = suggestVarName('result');
    assignin('caller', char(suggestedName), rand(3,3));
    fprintf('Created variable: %s\n', suggestedName);

    % Show updated workspace
    updatedVars = who;
    fprintf('\nUpdated workspace variables: %s\n', strjoin(updatedVars, ', '));
end

% Enhanced version with additional features
function newVarName = suggestVarNameAdvanced(baseName, varargin)
    % SUGGESTVARNAMEADVANCED Advanced variable name suggestion with options
    %
    % Usage:
    %   newVarName = suggestVarNameAdvanced('myVar')
    %   newVarName = suggestVarNameAdvanced('myVar', 'workspace', 'base')
    %   newVarName = suggestVarNameAdvanced('myVar', 'separator', '_')
    %   newVarName = suggestVarNameAdvanced('myVar', 'startNumber', 0)
    %   newVarName = suggestVarNameAdvanced('myVar', 'maxSuggestions', 5)
    %
    % Input:
    %   baseName - string or char array with the desired base variable name
    %
    % Optional Parameters:
    %   workspace      - 'base', 'caller', or 'auto' (default: 'auto')
    %   separator      - character to use between base name and number (default: '')
    %   startNumber    - starting number for suggestions (default: 1)
    %   maxSuggestions - maximum number of suggestions to try (default: 10000)
    %   returnAll      - return all available names up to maxSuggestions (default: false)
    %
    % Output:
    %   newVarName - string with suggested variable name(s)

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'baseName', @(x) ischar(x) || isstring(x));
    addParameter(p, 'workspace', 'auto', @(x) ismember(lower(x), {'base', 'caller', 'auto'}));
    addParameter(p, 'separator', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'startNumber', 1, @(x) isnumeric(x) && x >= 0);
    addParameter(p, 'maxSuggestions', 10000, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'returnAll', false, @islogical);
    parse(p, baseName, varargin{:});

    % Extract parameters
    baseName = string(p.Results.baseName);
    workspace = p.Results.workspace;
    separator = string(p.Results.separator);
    startNumber = p.Results.startNumber;
    maxSuggestions = p.Results.maxSuggestions;
    returnAll = p.Results.returnAll;

    % Determine which workspace to use
    if strcmp(workspace, 'auto')
        try
            evalin('caller', 'who');
            workspace = 'caller';
        catch
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
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:evalin:invalidWorkspace')
            existingVars = evalin('base', 'who');
        else
            rethrow(ME);
        end
    end

    % Convert existing variable names to strings for comparison
    existingVars = string(existingVars);

    % Check if base name exists
    if ~any(strcmp(existingVars, baseName))
        newVarName = baseName;
        if ~returnAll
            return;
        end
    end

    % Find available names
    availableNames = string.empty;
    counter = startNumber;

    while length(availableNames) < maxSuggestions
        if counter == 0 && ~any(strcmp(existingVars, baseName))
            % Special case: if startNumber is 0 and base name is available
            candidateName = baseName;
        else
            % Create candidate name with separator and number
            if counter == 0
                counter = 1;  % Skip 0 if base name exists
            end
            candidateName = baseName + separator + string(counter);
        end

        % Check if this candidate exists in the variable list
        if ~any(strcmp(existingVars, candidateName))
            availableNames(end+1) = candidateName;
            if ~returnAll
                newVarName = candidateName;
                return;
            end
        end

        counter = counter + 1;

        % Safety check to prevent infinite loops
        if counter > startNumber + maxSuggestions * 2
            break;
        end
    end

    if returnAll
        newVarName = availableNames;
    else
        if ~isempty(availableNames)
            newVarName = availableNames(1);
        else
            error('Could not find available variable name within the specified limits');
        end
    end
end

% Utility function to get all variable names in workspace
function varNames = getWorkspaceVars()
    % Returns a cell array of all variable names in caller's workspace
    varNames = evalin('caller', 'who');
end

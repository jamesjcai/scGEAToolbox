function i_resetrngseed
    % Create a simple GUI for selecting a random seed
    
    % Create the main figure window
    hFig = figure('Position', [500, 500, 300, 150], 'MenuBar', 'none', ...
                  'Name', 'Random Seed Selector', 'NumberTitle', 'off', ...
                  'Resize', 'off');
    
    % Create a label for the seed input
    uicontrol('Style', 'text', 'Position', [30, 80, 120, 25], ...
              'String', 'Enter seed value:');
    
    % Create an input field for the seed value
    hSeedInput = uicontrol('Style', 'edit', 'Position', [150, 80, 100, 25], ...
                           'BackgroundColor', 'white', 'String', '');
    
    % Create a button to set the seed
    uicontrol('Style', 'pushbutton', 'Position', [30, 30, 100, 30], ...
              'String', 'Set Seed', 'Callback', @setSeed);
    
    % Create a button to shuffle and set a random seed
    uicontrol('Style', 'pushbutton', 'Position', [150, 30, 100, 30], ...
              'String', 'Random Seed', 'Callback', @randomSeed);
    
    % Nested function to handle setting the seed
    function setSeed(~, ~)
        seedValue = str2double(get(hSeedInput, 'String'));
        if isnan(seedValue)
            errordlg('Please enter a valid numeric seed', 'Invalid Input');
        else
            rng(seedValue);
            msgbox(sprintf('Random seed set to: %d', seedValue), 'Success');
        end
    end

    % Nested function to handle setting a random seed
    function randomSeed(~, ~)
        rng('shuffle');
        seedValue = sum(100*clock);  % Generate a seed based on current time
        msgbox(sprintf('Random seed (shuffled) set to: %d', seedValue), 'Success');
    end
end

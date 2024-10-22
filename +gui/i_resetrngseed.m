function i_resetrngseed(~, ~)

    answer = questdlg("Set random seed.","", "Default Seed", "Random Seed", "Set Seed", "Default Seed");
    switch answer
        case 'Default Seed'
            rng("default");
            helpdlg(sprintf('Random seed set to default (%d).', 0), '');            
        case 'Random Seed'
            rng('shuffle');
            seedValue = sum(100*datevec(datetime));  % Generate a seed based on current time
            helpdlg(sprintf('Random seed (shuffled) set to: %d', seedValue), '');
        case 'Set Seed'
            seedValue = sum(100*datevec(datetime));  % Generate a seed based on current time
            seedValue = sprintf('%.0f', seedValue);
            seedValue = gui.i_inputnumk(seedValue, 0, 2^32, ...
                "Random number seed, specified as a nonnegative integer less than 2^32");
            if isempty(seedValue) || isnan(seedValue)
                errordlg('Please enter a valid numeric seed', 'Invalid Input');
            else
                rng(seedValue);
                helpdlg(sprintf('Random seed set to: %d', seedValue), '');
            end
    end

end
function i_resetrngseed(src, ~, needconfirm)

if nargin < 3, needconfirm = true; end

[parentfig, ~] = gui.gui_getfigsce(src);
    answer = gui.myQuestdlg(parentfig, "Set random seed.","", ...
        {'Default Seed', 'Random Seed', 'Set Seed'}, 'Default Seed');
    switch answer
        case 'Default Seed'
            rng("default");
            if needconfirm
                gui.myHelpdlg(parentfig, 'Random seed set to default.');
            end
        case 'Random Seed'
            rng('shuffle');
            seedValue = round(sum(100*datevec(datetime)));  % Generate a seed based on current time
            if needconfirm
                gui.myHelpdlg(parentfig, sprintf('Random seed (shuffled) set to: %d', seedValue));
            end
        case 'Set Seed'
            seedValue = sum(100*datevec(datetime));  % Generate a seed based on current time
            seedValue = sprintf('%.0f', seedValue);
            seedValue = gui.i_inputnumk(seedValue, 0, 2^32, ...
                "Random number seed, specified as a " + ...
                "nonnegative integer less than 2^32", parentfig);
            if isempty(seedValue) || isnan(seedValue)
                gui.myErrordlg(parentfig, ['Please enter a valid ' ...
                    'numeric seed'], 'Invalid Input');
            else
                rng(seedValue);
                if needconfirm
                    gui.myHelpdlg(parentfig, ...
                        sprintf('Random seed set to: %d', ...
                        seedValue));
                end
            end
    end

end
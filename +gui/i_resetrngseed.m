function ok = i_resetrngseed(src, ~, needconfirm)
%I_RESETRNGSEED Ask for a random seed and set it.
%
%   gui.i_resetrngseed(src) offers the default seed, a shuffled one, or a
%   seed the user types, sets it with RNG, and says which seed is now in
%   use. NEEDCONFIRM = false skips that last message.
%
%   OK is false when the dialog was cancelled or closed, so a caller can
%   stop instead of running on with the generator unchanged.
%
%   The seed reported after a shuffle is the one RNG chose, read back from
%   it, so a run can be repeated from the number shown.

if nargin < 3, needconfirm = true; end
ok = false;
maxSeed = 2^32 - 1;   % RNG takes a nonnegative integer below 2^32

[parentfig, ~] = gui.gui_getfigsce(src);
answer = gui.myQuestdlg(parentfig, "Set random seed.", "", ...
    {'Default Seed', 'Random Seed', 'Set Seed'}, 'Default Seed');
switch answer
    case 'Default Seed'
        rng("default");
        msg = 'Random seed set to default.';
    case 'Random Seed'
        rng('shuffle');
        s = rng;
        msg = sprintf('Random seed (shuffled) set to: %d', s.Seed);
    case 'Set Seed'
        s = rng;
        seedValue = gui.i_inputnumk(s.Seed, 0, maxSeed, ...
            "Random number seed, specified as a " + ...
            "nonnegative integer less than 2^32", parentfig);
        % Empty on Cancel, and on a bad value I_INPUTNUMK has already
        % said what was wrong.
        if isempty(seedValue), return; end
        rng(seedValue);
        msg = sprintf('Random seed set to: %d', seedValue);
    otherwise
        % Dismissed: leave the generator as it was.
        return;
end
ok = true;
if needconfirm
    gui.myHelpdlg(parentfig, msg);
end
end

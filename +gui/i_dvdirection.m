function [direction, labels] = i_dvdirection(parentfig, direction)
%I_DVDIRECTION Ask how DV genes are split into up and down, and name them.
%
%   [DIRECTION, LABELS] = GUI.I_DVDIRECTION(PARENTFIG) asks the user and
%   returns 'mean' or 'deviation' (see SC_DVG), or '' if the dialog was
%   dismissed.
%
%   [DIRECTION, LABELS] = GUI.I_DVDIRECTION([], DIRECTION) returns the
%   labels for DIRECTION without asking.
%
%   LABELS has fields:
%       Up, Dn       - list names, e.g. 'Up-regulated'
%       UpTag, DnTag - short sheet-name prefixes, e.g. 'Up'
%       FileTag      - file-name suffix: '' for 'mean', '_devsign' otherwise,
%                      so the two kinds of run do not overwrite each other

if nargin < 2
    % The explanation goes in the message, not on the buttons: UICONFIRM
    % sizes each button to its label, so a sentence-long option makes a
    % button wide enough to push the dialog off the screen.
    optMean = 'Mean expression';
    optDev = 'Deviation from curve';
    message = {'How should DV genes be split into up and down?'; ''; ...
        [optMean ' -- up = higher mean in group 1']; ...
        [optDev ' -- up = more variable in group 1']};
    answer = gui.myQuestdlg(parentfig, message, ...
        'DV Direction', {optMean, optDev}, optMean);
    switch answer
        case optMean
            direction = 'mean';
        case optDev
            direction = 'deviation';
        otherwise
            % Dismissed.
            direction = '';
    end
end

switch direction
    case 'mean'
        labels = struct('Up', 'Up-regulated', 'Dn', 'Down-regulated', ...
            'UpTag', 'Up', 'DnTag', 'Dn', 'FileTag', '');
    case 'deviation'
        labels = struct('Up', 'Variability increasing', ...
            'Dn', 'Variability decreasing', ...
            'UpTag', 'VInc', 'DnTag', 'VDec', 'FileTag', '_devsign');
    otherwise
        % Dismissed, or not a direction SC_DVG knows.
        labels = struct([]);
end
end

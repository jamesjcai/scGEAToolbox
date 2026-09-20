function [tf] = i_confirmonce(parentfig, key, message, title)
%I_CONFIRMONCE Ask for confirmation once per session, then take it as given.
%
%   tf = gui.i_confirmonce(parentfig, key, message)
%   tf = gui.i_confirmonce(parentfig, key, message, title)
%   gui.i_confirmonce([])                                  ask again
%
%   Shows MESSAGE with Continue/Cancel and returns true when the user
%   continues. KEY names the notice; once it has been confirmed, every later
%   call with the same KEY returns true without a dialog.
%
%   This is for a notice the user has to read once to work safely and never
%   again - "these labels are not saved" before a one-time analysis they may
%   run on twenty selections in a row. Cancelling is not remembered: the next
%   call asks again, because the answer was no.
%
%   Do not use it for a question whose answer can differ between calls, or for
%   one whose consequence is worth re-reading, such as overwriting data. Those
%   belong in GUI.MYQUESTDLG.
%
%   The record lives in a persistent variable, so it lasts as long as the
%   MATLAB session and is shared by every figure. CLEAR FUNCTIONS resets it, as
%   does gui.i_confirmonce([]).
%
%   See also gui.myQuestdlg, gui.i_sessionmarkers.

persistent confirmed

if isempty(confirmed), confirmed = strings(0, 1); end
if nargin < 2
    confirmed = strings(0, 1);
    tf = false;
    return;
end
if nargin < 4, title = ''; end

key = string(key);
if any(confirmed == key)
    tf = true;
    return;
end

answer = gui.myQuestdlg(parentfig, message, title, ...
    {'Continue', 'Cancel'}, 'Continue');
tf = strcmp(answer, 'Continue');

% Only a yes is remembered. A no means the user has not agreed to anything,
% so the next call has to ask.
if tf, confirmed(end+1, 1) = key; end
end

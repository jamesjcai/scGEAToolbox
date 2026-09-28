function [answer] = i_selvariabletype(y, parentfig)
if nargin<2, parentfig = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

[c] = findgroups(y);
n = max(c);
if n < 20
    deft = 'Categorical/Discrete';
else
    deft = 'Numerical/Continuous';
end
answer = gui.myQuestdlg(parentfig, 'What is the variable type?', '', ...
{'Categorical/Discrete', ...
'Numerical/Continuous', 'Unknown'}, deft);
if isempty(answer)
    % Closed: a cancel. It used to become 'Unknown', which the app treats
    % as numerical, so a dismissed dialog went on to recolor the cells.
    answer = '';
    return;
end


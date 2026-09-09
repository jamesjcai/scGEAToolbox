function k = i_inputnumk(defaultk, a, b, descstr, parentfig)

if nargin < 5, parentfig = []; end
if nargin < 4 || isempty(descstr), descstr = 'Enter a number'; end
if nargin < 3 || isempty(b), b = 100; end
if nargin < 2 || isempty(a), a = 1; end
if nargin < 1, defaultk = 10; end
if ~isempty(parentfig)
    gui.i_raisefig(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
k = [];
prompt = {sprintf('%s (%d..%d):', ...
descstr, a, b)};
dlgtitle = '';
dims = [1, 50];
if isnumeric(defaultk)
    definput = {sprintf('%d', defaultk)};
else
    definput = {defaultk};
end

if gui.i_isuifig(parentfig)
    % answer = gui.myInputdlg(prompt, dlgtitle, definput, parentfig);
    answer = gui.i_inputdlg(prompt, definput, parentfig);
else
    answer = inputdlg(prompt, dlgtitle, dims, definput);
end
if isempty(answer), return; end
k = round(str2double(cell2mat(answer)));
if isnan(k) || k < a || k > b
    k = [];
    % Route through myErrordlg so a uifigure parent gets an in-figure uialert
    % instead of a separate modal window that steals the window ordering.
    gui.myErrordlg(parentfig, ...
        sprintf('Enter a whole number between %d and %d.', a, b), '', false);
    return;
end
end

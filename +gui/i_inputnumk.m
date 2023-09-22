function k = i_inputnumk(defaultk, a, b, descstr)
if nargin < 4, descstr = 'a number'; end
if nargin < 3, b = 100; end
if nargin < 2, a = 1; end
if nargin < 1, defaultk = 10; end
k = [];
prompt = {sprintf('Enter %s (%d..%d):', ...
    descstr, a, b)};
dlgtitle = '';
dims = [1, 45];
definput = {sprintf('%d', defaultk)};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer), return; end
k = round(str2double(cell2mat(answer)));
if isnan(k) || k < a || k > b
    k = [];
    uiwait(errordlg('Invalid number.'));
    return;
end
end

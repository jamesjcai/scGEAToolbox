function k = i_inputnumk(defaultk)
    if nargin<1, defaultk=10; end
    k=[];
    prompt = {'Enter number of k = (2..100):'};
    dlgtitle = 'Input k';
    dims = [1 45];
    definput = {sprintf('%d',defaultk)};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(answer), return; end
    k = round(str2double(cell2mat(answer)));
    if isnan(k) || k < 2 || k > 100
        k=[];
        uiwait(errordlg('Invalid K'));
        return;
    end
end

function k = i_inputnumk
    k=[];
    prompt = {'Enter number of clusters K=(2..100):'};
    dlgtitle = 'Input K';
    dims = [1 45];
    definput = {'10'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if isempty(answer), return; end
    k = round(str2double(cell2mat(answer)));
    if isnan(k) || k < 2 || k > 100
        k=[];
        uiwait(errordlg('Invalid K'));
        return;
    end
end

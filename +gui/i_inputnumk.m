    function k = i_inputnumk
        prompt = {'Enter number of clusters K=(2..50):'};
        dlgtitle = 'Input K';
        dims = [1 45];
        definput = {'10'};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer)
            k=[];
            return
        end
        k = round(str2double(cell2mat(answer)));
    end

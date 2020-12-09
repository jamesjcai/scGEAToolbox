function [numfig]=gui_inputdlg

    prompt = {'Enter number of genes (1-50)'};
    dlgtitle = ' ';
    answer = inputdlg(prompt,dlgtitle,[1 50],{'10'});
    if isempty(answer)
        return; 
    else
        try
            numfig=str2double(answer{1});
        catch ME
            errordlg(ME.message);
            return;
        end
    end
    if ~(numfig>0 && numfig<=50)
        errordlg('Invalid number of figures');
        return;
    end

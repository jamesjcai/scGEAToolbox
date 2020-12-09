function [numfig]=gui_inputdlg

    prompt = {'Enter number of panels (1-10)'};
    dlgtitle = 'Panel of 9 Genes';
    answer = inputdlg(prompt,dlgtitle,[1 40],{'1'});
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
    if ~(numfig>0 && numfig<=10)
        errordlg('Invalid number of figures');
        return;
    end

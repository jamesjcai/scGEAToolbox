function [glist]=gui_inputgenelist    
    prompt = {'Paste Gene List:'};
    dlgtitle = 'Input Genes';    
    
[answer] = inputdlg(prompt,...
    dlgtitle,[20 40],...
    {sprintf('Gene1\nGene2')});
if isempty(answer)
    glist=[];
    return;
end
if iscell(answer)
    answer=string(answer{1});
end
try
    glist=answer;
catch ME
    errordlg(ME.message);
    return;
end
end


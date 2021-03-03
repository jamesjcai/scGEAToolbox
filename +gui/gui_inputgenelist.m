function [glist]=gui_inputgenelist(glist0)
    if nargin<1
        glist0=["Gene1";"Gene2";"Gene3";...
            "Gene4";"Gene5";"Gene6";...
            "Gene7";"Gene8";"Gene9";"Gene10"];
    end
    s=sprintf('%s\n',glist0);
    s=s(1:end-1);
    prompt = {'Paste Gene List:'};
    dlgtitle = 'Input Genes';    
[answer] = inputdlg(prompt,...
    dlgtitle,[20 40],{s});
glist=[];
if isempty(answer), return; end
if iscell(answer)
try
    glist=string(answer{1});
    glist=strrep(glist,'"','');
    glist=strip(glist,'both','"');
    glist=strtrim(glist);
catch ME
    errordlg(ME.message);
    return;
end
end
end


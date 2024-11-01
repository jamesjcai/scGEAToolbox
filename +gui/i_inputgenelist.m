function [glist] = i_inputgenelist(glist0, allowspace)
if nargin < 2, allowspace = false; end
if nargin < 1
    glist0 = ["Gene1"; "Gene2"; "Gene3"; ...
        "Gene4"; "Gene5"; "Gene6"; ...
        "Gene7"; "Gene8"; "Gene9"; "Gene10"];
end
s = sprintf('%s\n', glist0);
s = s(1:end-1);
prompt = {'Paste List:'};
dlgtitle = 'Input List';
[answer] = inputdlg(prompt, dlgtitle, [15, 80], {s});
glist = [];
if isempty(answer), return; end
if iscell(answer)
    try
        glist = string(answer{1});
        glist = strrep(glist, '"', '');
        glist = strrep(glist, '''', '');
        glist = strip(glist, 'both', '"');
        glist = strtrim(glist);

        if allowspace
            spset = {',', ';'};
        else
            spset = {' ', ',', ';'};
        end
        i = contains(glist, spset);
        if any(i)
            g1 = glist(~i);
            g2 = glist(i);
            g3 = [];
            for k = 1:length(g2)
                gx = strsplit(g2(k), spset);
                g3 = [g3; gx.'];
            end
            glist = [g1; g3];
            glist = strtrim(glist);
        end
    catch ME
        errordlg(ME.message);
        return;
    end
end
end

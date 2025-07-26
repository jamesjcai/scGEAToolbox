function [markerg] = i_get_panglaodbmarkers
    markerg = "";
    try
        pw1 = fileparts(mfilename('fullpath'));        
        pth = fullfile(pw1, '..', 'external', 'fun_alona_panglaodb','marker_hs.mat');
        load(pth,'Tw');
        markerg1 = unique(string(Tw.Var1));
        % idx(ismember(upper(obj.g), upper(markerg))) = true;
        pth = fullfile(pw1, '..', 'external', 'fun_alona_panglaodb','marker_mm.mat');
        load(pth,'Tw');
        markerg2 = unique(string(Tw.Var1));
        % idx(ismember(upper(obj.g), upper(markerg))) = true;
        markerg = unique(upper([markerg1; markerg2]));
        %fprintf('EMBEDCELLS: %d additional marker genes included.\n', sum(idx(numhvg+1:end)));
    catch
        %warning('EMBEDCELLS: marker genes not included.');
    end
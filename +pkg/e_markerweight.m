function [Tw] = e_markerweight(Tm)

    s = upper(string(Tm.PositiveMarkers));
    S = [];
    for k = 1:length(s)
        sk = s(k);
        a = strsplit(sk, ',');
        a = strtrim(a);
        if strlength(a(end)) == 0 || isempty(a(end))
            a = a(1:end-1);
        end        
        S = [S, a];
    end
    
    %%
    N = length(S);
    t = tabulate(S);
    f = cell2mat(t(:, 3));
    if max(f) - min(f) < eps
        w = ones(N, 1);
    else
        w = 1 + sqrt((max(f) - f)/(max(f) - min(f)));
    end
    genelist = string(t(:, 1));
    Tw = table(genelist, w);
    Tw.Properties.VariableNames = {'Var1', 'Var2'};


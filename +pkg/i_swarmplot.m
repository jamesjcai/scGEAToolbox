function i_swarmplot(d, c)

if ~isstring(c), c = string(c); end

[c, cL] = grp2idx(c);
swarmchart(c, d);
end

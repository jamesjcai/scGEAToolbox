function [new]=i_makeuniquename(orig)
new = orig;
uniqStr = unique(orig);
for ii = 1:length(uniqStr)
    %idx = strmatch(uniqStr{ii}, orig, 'exact');
    idx = find(strcmp(uniqStr{ii}, orig));
    for jj = 2:length(idx)
        new{idx(jj)} = [orig{idx(jj)}, '_', num2str(jj)];
    end
end
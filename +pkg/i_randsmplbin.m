function i=i_randsmplbin(idx,bin)
    i=[];
    binpool=bin(idx);
    ubinpool=unique(binpool);
    for k=1:length(ubinpool)
        tmpi=binpool==ubinpool(k);
        i=[i; randsample(find(bin==ubinpool(k)),sum(tmpi))];
    end
end
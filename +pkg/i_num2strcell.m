function g=i_num2strcell(d,prefix)
    if isscalar(d), d=1:d; end
    sf=sprintf('%s%%d,',prefix);
    g=textscan(sprintf(sf,d),'%s','delimiter',',');
    g=string(g{1});
end
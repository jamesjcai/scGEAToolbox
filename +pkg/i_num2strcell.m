function [gs,gc]=i_num2strcell(d,prefix)
    if nargin<2, prefix=''; end
    if isscalar(d), d=1:d; end
    sf=sprintf('%s%%d,',prefix);
    gc=textscan(sprintf(sf,d),'%s','delimiter',',');
    gc=gc{1};
    gs=string(gc);
end
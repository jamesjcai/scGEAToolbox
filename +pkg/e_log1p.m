function [x] = e_log1p(x)


if ~issparse(x)
    x = log(1+x);
else
    fhandle = @(x) x + 1;
    x = spfun(fhandle, x);
    x = spfun(@log, x);
end

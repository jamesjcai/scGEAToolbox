function x = getUnifPoints(df,n,h1,a,b,h2)
%getUnifPoints returns evenly-spaced points according to given speed/density
%function (df, which is a function handle)
%   n is number of evenly-spaced points desired
%   h1 is desired level of (horizontal) discretization of domain [a,b]
%   h2 is desired level of (vertical) discretization of density values

m = ceil((b-a)/h1)+1;
h1 = (b-a)/(m-1);
bins = [];
min_val = df(a);

for i=0:m-1
    val = df(a+i*h1);
    min_val = min(val,min_val);
    num_bins = floor(val/h2);
    bins = [bins,(a+i*h1)*ones(1,num_bins)];
end
tot_bins = length(bins);
del = floor((tot_bins-1)/(n-1));
x = bins(1:del:end)';

end


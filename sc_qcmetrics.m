function [nmad] = sc_qcmetrics(X)
%QC metrics

%{
https://alona.panglaodb.se/add.html

nmad - this field gives the number of median absolute deviations (MADs) 
below the median of two quality metrics for which a cell is removed.
The two quality metrics are the log of the library size and 
the log of the number of expressed genes.
It is enough for a cell to be below one of the metrics for it to be removed.
The default value is 3.
Give a higher value to remove less cells and a lower value to remove more cells

pct - Use genes expressed in at least X% of the cells
%}
nc = size(X, 2);
pct = sum(X > 0, 2) ./ nc; % 1%
reads_per_cell = sum(X, 2)'; % reads per cell
nmad = mad(log(sum(X, 2))-log(sum(X > 0, 2)), 1);

end
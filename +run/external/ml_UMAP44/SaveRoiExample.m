function SaveRoiExample(key, roi, name, reduction, args)
rows=RoiUtil.GetRows(roi, reduction);
[R,nReducedCols]=size(reduction);
fprintf(['roi ="%s"/%s, reduction size=%d%d, ' ...
    '# selected rows=%d,\n\tmedians ' ...
    'on region''s unreduced data ='], name, ...
    key, R, nReducedCols, sum(rows));
subset=args.unreduced_data(rows, :);
mdns=median(subset,1);
cols=args.unreduced_column_names;
nUnreducedCols=length(cols);
for c=1:nUnreducedCols
    fprintf('\t%s=%s\n', cols{c}, num2str(mdns(c)));
end
end

LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};
Age = [38;43;38;40;49];
Smoker = logical([1;0;1;0;1]);
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
T = table(LastName,Age,Smoker,Height,Weight);

fname = "aaa";
parquetwrite(fname, T);
tbl = parquetread(fname);

% https://www.linkedin.com/posts/the-mathworks_2_using-matlab-with-python-activity-7363510088839434243-4haI

% from python
% df = pandas.read_parquet(fname)
% pandas.Dataframe.to_parquet(df)
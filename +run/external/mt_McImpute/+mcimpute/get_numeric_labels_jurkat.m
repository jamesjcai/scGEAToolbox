function [ actual_labels ] = get_numeric_labels_jurkat(dataname, data_dir  )
%GET_NUMERIC_LABELS_JURKAT Summary of this function goes here
%   Detailed explanation goes here

fid=fopen([data_dir 'annotations/' dataname '_annotations.csv']);
gt=textscan(fid,'%s','Delimiter',',');
gt_labels=gt{1};
keySet =   { 'species1','species2'};
valueSet = [1,2];
mapi=containers.Map(keySet,valueSet); %access labels using map(gt{i})
actual_labels=[];
for i=1:length(gt_labels)
actual_labels=[actual_labels;mapi(gt_labels{i})];
end
end


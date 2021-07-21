function [labels] = get_numeric_labels_Zeisel(dataname, data_dir  )

load([data_dir 'annotations/' dataname '_annotations.mat']);
gt_labels=anno;
save('Temp/Zeisel_labelNames.mat','gt_labels');
keySet =     {'interneurons','s1pyramidal', 'ca1pyramidal', 'oligodendrocytes','microglia','endothelial' , 'astrocytes','ependymal' ,'mural' };
valueSet = [1,2,3,4,5,6,7,8,9];
mapi=containers.Map(keySet,valueSet); %access labels using map(gt{i})
labels=[];
for i=1:length(gt_labels)
labels=[labels;mapi(gt_labels{i})];
end

end


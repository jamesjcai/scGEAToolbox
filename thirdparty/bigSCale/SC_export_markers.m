%SC_export_markers
delete('./../data/output_export_markers/markers_hard.xlsx')
delete('./../data/output_export_markers/markers_soft.xlsx')
copyfile('./../data/output_export_markers/template.xlsx','./../data/output_export_markers/markers_hard.xlsx')
copyfile('./../data/output_export_markers/template.xlsx','./../data/output_export_markers/markers_soft.xlsx')

%info=kgXref_human_final;
%info=repmat(gene_names,1,3);
info=all_ensambl_genes;

for k=1:length(markers_hard)

    A=info(markers_hard{k},:);
    if isrow(scores_hard{k})
        scores_hard{k}=scores_hard{k}';
    end
    B=mio_mat2cell(scores_hard{k});
    C=mio_mat2cell(FC_hard{k});
    if isrow(FC_hard{k})
        C=C';
    end
    %header={'Gene ID' 'Gene name' 'Synonims' 'Gene type' 'Description' 'chr' 'Marker score'};
    header={'Gene name' 'Gene type' 'Description' 'Marker score' 'Marker log2(FC)'};
    if length(scores_hard{k})>0
    xlswrite('./../data/output_export_markers/markers_hard.xlsx',vertcat(header,[A B C]),sprintf('C%g',k),'A1');
    end
end

for k=1:length(markers_soft)
    
    A=info(markers_soft{k},:);
    if isrow(scores_soft{k})
        scores_soft{k}=scores_soft{k}';
    end
    B=mio_mat2cell(scores_soft{k});
    C=mio_mat2cell(FC_soft{k});
        if isrow(FC_soft{k})
        C=C';
    end
    %header={'Gene ID' 'Gene name' 'Synonims' 'Gene type' 'Description' 'chr' 'Marker score'};
    header={'Gene name' 'Gene type' 'Description' 'Marker score' 'Marker log2(FC)'};
    if length(scores_soft{k})>0
    xlswrite('./../data/output_export_markers/markers_soft.xlsx',vertcat(header,[A B C]),sprintf('C%g',k),'A1');
    end
end
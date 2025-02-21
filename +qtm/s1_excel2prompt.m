t = readtable('test.xlsx','Sheet','Up_250_GO_BP');
T = t(:,[3 7]);
T.Properties.VariableNames={'Function Term','Genes'};
T.Genes = strrep(T.Genes,',', ', ');
s = jsonencode(table2struct(T));

%{
t = readtable('DEVP_DP_C2_x08w_vs_x14w_AlphaCells.xlsx','Sheet','Up-regulated');
t = readtable('DEVP_DP_C2_x08w_vs_x14w_AlphaCells.xlsx','Sheet','Down-regulated');
T = t(:, 1);
s = jsonencode(table2struct(T));

t = readtable('DEVP_DV_x08w_vs_x14w_AlphaCells.xlsx','Sheet','Up_250_GO_BP');
T = t(:,[1 3 6 7 8]);
s = jsonencode(table2struct(T))
%}
% s='[{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Negative Regulation Of Cell Cycle (GO:0045786)","CombinedScore":30.002437542843889,"OverlappingGenes":"RNF167,CDKN1A,NUPR1,IPO5,RHOB","AdjustedP_value":0.52379660316216059},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Mitochondrial Gene Expression (GO:0140053)","CombinedScore":12.321933466777262,"OverlappingGenes":"MRPS12,MRPL38,FOXO3,MRPL13,MRPL43,MRPL54,DAP3,MRPL33","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Regulation Of Cell Cycle (GO:0051726)","CombinedScore":9.2576864622077331,"OverlappingGenes":"SDCBP,RNF167,CDKN1A,SON,GADD45A,CDC26,MORF4L2,STAT3,NUPR1,EIF4G2,RHOB,GADD45G","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Regulation Of Neuron Apoptotic Process (GO:0043523)","CombinedScore":12.802762668081101,"OverlappingGenes":"EGLN2,NUPR1,APOE,IL6ST,FOXO3","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Mitochondrial Translation (GO:0032543)","CombinedScore":8.81937831079096,"OverlappingGenes":"MRPS12,MRPL38,MRPL13,MRPL43,MRPL54,DAP3,MRPL33","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Negative Regulation Of Cell Population Proliferation (GO:0008285)","CombinedScore":7.0152985541582469,"OverlappingGenes":"KLF10,PBRM1,CDKN1A,IFITM1,RBBP4,LMNA,IGFBP7,APOE,NUPR1,B2M,EMD","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Cellular Response To Starvation (GO:0009267)","CombinedScore":9.5257289543113135,"OverlappingGenes":"SH3GLB1,KLF10,GABARAPL1,CDKN1A,MAP1LC3A,FOXO3","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Regulation Of Cell Migration (GO:0030334)","CombinedScore":6.7456362041280347,"OverlappingGenes":"SPAG9,SDCBP,CLDN4,IFITM1,SERPINE2,RBBP4,TPM1,LMNA,STAT3,FOXO3,RHOB","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Response To Cytokine (GO:0034097)","CombinedScore":9.080444376009547,"OverlappingGenes":"IFITM3,IFITM1,TIMP3,IL6ST,DDOST","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Positive Regulation Of Protein Catabolic Process (GO:0045732)","CombinedScore":9.080444376009547,"OverlappingGenes":"NDUFA13,EGLN2,FMR1,NUPR1,APOE","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Central Nervous System Development (GO:0007417)","CombinedScore":6.808364982248917,"OverlappingGenes":"NR4A2,DHX30,GIT2,RBBP4,UTP3,STAT3,RAB18,THOC6","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Regulation Of mRNA Splicing, Via Spliceosome (GO:0048024)","CombinedScore":7.3360449932494713,"OverlappingGenes":"KHDRBS3,SON,FMR1,CELF3,CIRBP,HNRNPA0","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Negative Regulation Of Apoptotic Process (GO:0043066)","CombinedScore":4.4486346314721228,"OverlappingGenes":"ARL6IP1,HSPB1,THOC6,DNAJA1,NR4A2,SOCS3,KRT18,SON,TMBIM4,APOE,NUPR1,NAA15,IL6ST","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Axonogenesis (GO:0007409)","CombinedScore":6.3783335037666458,"OverlappingGenes":"ENAH,USP9X,NRXN1,S100A6,PAK3","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Negative Regulation Of Gene Expression (GO:0010629)","CombinedScore":4.64235550899219,"OverlappingGenes":"SRGN,OLFM1,FMR1,AGO2,STAT3,CNPY2,APOE,ATP2B1,PAIP2,MRPL13","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Negative Regulation Of Cell Motility (GO:2000146)","CombinedScore":5.868926747849164,"OverlappingGenes":"IFITM1,RBBP4,TPM1,FOXO3,RHOB","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Regulation Of Apoptotic Process (GO:0042981)","CombinedScore":3.6959861029919869,"OverlappingGenes":"EGLN2,ARL6IP1,GADD45A,HSPB1,FOXO3,RHOB,GADD45G,THOC6,DNAJA1,SOCS3,OLFM1,KRT18,SON,TMBIM4,MORF4L2,NUPR1,NAA15,IL6ST","AdjustedP_value":0.55500058821099107},{"GeneSetLibrary":"GO_Biological_Process_2023","TermName":"Cytokine-Mediated Signaling Pathway (GO:0019221)","CombinedScore":4.9938634421547254,"OverlappingGenes":"PLVAP,STAT3,P4HB,IL6ST,FOXO3","AdjustedP_value":0.55500058821099107}]';

chat = ollamaChat("deepseek-r1", TimeOut = 1200);
prompt1 = "enrichr is a gene function enrichment analysis service. I will give you an output of enrichr analysis below. Please summarize the results in text. The results mixed enriched terms of biological processes and molecular functions. Please write the summary in paragraph(s). Do not use bullet points. ";
prompt2 = "Here is the output of enrichr: " + s;
feedbk = generate(chat, prompt1 + prompt2);        

import mlreportgen.dom.*
doc = Document('test', 'docx');
open(doc);
para = Paragraph(feedbk);
append(doc, para);

close(doc);
rptview('test', 'docx');


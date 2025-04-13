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


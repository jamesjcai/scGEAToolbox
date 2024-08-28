function [response] = mt_Enrichr(genelist, backgroundlist, genesets, minumgene)


if nargin < 4, minumgene = 5; end

if nargin < 3
    genesets = ["GO_Biological_Process_2023",...
                "GO_Molecular_Function_2023"]; 
end

if nargin < 2
    backgroundlist = ["NSUN3","POLRMT","NLRX1","SFXN5","ZC3H12C", ...
    "SLC25A39","ARSG",...
	"DEFB29","PCMTD2","ACAA1A","LRRC1","2810432D09RIK","SEPHS2",...
	"SAC3D1","TMLHE","LOC623451","TSR2","PLEKHA7","GYS2","ARHGEF12",...
	"HIBCH","LYRM2","ZBTB44","ENTPD5","RAB11FIP2","LIPT1",...
	"INTU","ANXA13","KLF12","SAT2","GAL3ST2","VAMP8","FKBPL",...
	"AQP11","TRAP1","PMPCB","TM7SF3","RBM39","BRI3","KDR","ZFP748",...
	"NAP1L1","DHRS1","LRRC56","WDR20A","STXBP2","KLF1","UFC1",...
	"CCDC16","9230114K14RIK","RWDD3","2610528K11RIK","ACO1",...
	"CABLES1", "LOC100047214","YARS2","LYPLA1","KALRN","GYK",...
	"ZFP787","ZFP655","RABEPK","ZFP650","4732466D17RIK","EXOSC4",...
	"WDR42A","GPHN","2610528J11RIK","1110003E01RIK","MDH1","1200014M14RIK",...
	"AW209491","MUT","1700123L14RIK","2610036D13RIK",...
	"PHF14", "RBM3", "MSL1", "PHF21A", "ARL10", "INSR", "JADE2", ...
	"P2RX7", "LINC00662", "CCDC101", "PPM1B", "KANSL1L", "CRYZL1", ...
	"ANAPC16", "TMCC1","CDH8", "RBM11", "CNPY2", "HSPA1L", "CUL2", ...
	"PLBD2", "LARP7", "TECPR2", "ZNF302", "CUX1", "MOB2", "CYTH2", ...    
	"SEC22C", "EIF4E3", "ROBO2", "ADAMTS9-AS2", "CXXC1", "LINC01314", "ATF7", ...
	"ATP5F1""COX15","TMEM30A","NSMCE4A","TM2D2","RHBDD3","ATXN2","NFS1",...
	"3110001I20RIK","BC038156","C330002I19RIK","ZFYVE20","POLI","TOMM70A",...
	"LOC100047782","2410012H22RIK","RILP","A230062G08RIK",...
	"PTTG1IP","RAB1","AFAP1L1", "LYRM5","2310026E23RIK",...
	"SLC7A6OS","MAT2B","4932438A13RIK","LRRC8A","SMO","NUPL2"];
    backgroundlist = [];
end
if nargin < 1
    genelist = ["PHF14", "RBM3", "MSL1", "PHF21A", "ARL10", "INSR", "JADE2", "P2RX7", ...
         "LINC00662", "CCDC101", "PPM1B", "KANSL1L", "CRYZL1", "ANAPC16", "TMCC1", ...
         "CDH8", "RBM11", "CNPY2", "HSPA1L", "CUL2", "PLBD2", "LARP7", "TECPR2", ...
         "ZNF302", "CUX1", "MOB2", "CYTH2", "SEC22C", "EIF4E3", "ROBO2", ...
         "ADAMTS9-AS2", "CXXC1", "LINC01314", "ATF7", "ATP5F1"];
end

import matlab.net.URI
import matlab.net.http.RequestMessage
import matlab.net.http.io.MultipartProvider
import matlab.net.http.io.MultipartFormProvider


base_url = 'https://maayanlab.cloud/Enrichr/addList';
formData = MultipartFormProvider('list', strjoin(genelist, newline), ...
           'description', 'genelist');
request = RequestMessage('post', [], formData);
response = send(request, URI(base_url));
res = jsondecode(response.Body.Data);
user_list_id = res.userListId;


gene_set_library = genesets(1); % "KEGG_2015";
ENRICHR_URL = "https://maayanlab.cloud/Enrichr/enrich";
query_string = sprintf("?userListId=%d&backgroundType=%s", ...
               user_list_id, gene_set_library); 
url = ENRICHR_URL + query_string;
response = jsondecode(convertCharsToStrings(char(webread(url))));


res = response.(gene_set_library);
isok = false(length(res),1);
for k = 1:length(res)
    if size(res{k}{6},1) >= minumgene
        isok(k) = true;
    end
end
response = res(isok);



%{
ENRICHR_URL = "https://maayanlab.cloud/Enrichr/export";
query_string = sprintf("?userListId=%d&filename=%s&backgroundType=%s", ...
    user_list_id, 'example_enrichment', gene_set_library);
url = ENRICHR_URL + query_string;
response = webread(url); % jsondecode(convertCharsToStrings(char(webread(url))));
response = strsplit(convertCharsToStrings(char(response)),'\n');
%}



%{
import matlab.net.URI
import matlab.net.http.RequestMessage
import matlab.net.http.io.MultipartProvider
import matlab.net.http.io.MultipartFormProvider

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList';
genes = {'PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2', 'P2RX7', ...
         'LINC00662', 'CCDC101', 'PPM1B', 'KANSL1L', 'CRYZL1', 'ANAPC16', 'TMCC1', ...
         'CDH8', 'RBM11', 'CNPY2', 'HSPA1L', 'CUL2', 'PLBD2', 'LARP7', 'TECPR2', ...
         'ZNF302', 'CUX1', 'MOB2', 'CYTH2', 'SEC22C', 'EIF4E3', 'ROBO2', ...
         'ADAMTS9-AS2', 'CXXC1', 'LINC01314', 'ATF7', 'ATP5F1'};
genes_str = strjoin(genes, newline);
description = 'Example gene list';

formData = MultipartFormProvider('list', genes_str, 'description', description);
request = RequestMessage('post', [], formData);
response = send(request, URI(ENRICHR_URL));

res = jsondecode(response.Body.Data);
user_list_id = res.userListId;


gene_set_library = "KEGG_2015";
ENRICHR_URL = "https://maayanlab.cloud/Enrichr/enrich";
query_string = sprintf("?userListId=%d&backgroundType=%s", ...
               user_list_id, gene_set_library); 
res = jsondecode(convertCharsToStrings(char(webread(ENRICHR_URL + query_string))));


% https://www.mathworks.com/matlabcentral/answers/2148269-try-to-call-the-rest-apis-provided-by-enrichr-from-matlab-but-webwrite-does-not-work#answer_1506069

%} 


%{
base_url = "https://maayanlab.cloud/speedrichr/api/addList";
formData = MultipartFormProvider('list', strjoin(genelist, newline), ...
           'description', 'genelist');
request = RequestMessage('post', [], formData);
response = send(request, URI(base_url));
res = response.Body.Data;
user_list_id = res.userListId;

base_url = "https://maayanlab.cloud/speedrichr/api/addbackground";
formData = MultipartFormProvider('background', strjoin(backgroundlist, newline));
request = RequestMessage('post', [], formData);
response = send(request, URI(base_url));
background_id = response.Body.Data.backgroundid;

base_url = "https://maayanlab.cloud/speedrichr/api/backgroundenrich";

formData = MultipartFormProvider('userListID', num2str(user_list_id), ...
            'backgroundid', background_id, ...
            'backgroundType', 'ChEA_2022');

request = RequestMessage('post', [], formData);
response = send(request, URI(base_url));
%}





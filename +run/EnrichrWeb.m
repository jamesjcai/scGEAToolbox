base_address = "https://maayanlab.cloud/Enrichr/";
sites_base_address = "https://maayanlab.cloud/";
sites = ["Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"];

url = strcat(sites_base_address, sites(1), "/", "datasetStatistics");
data=webread(url);
data.statistics{10}.libraryName

%%
url = strcat(base_address,"enrich");
uri = matlab.net.URI(url);
res = webwrite(uri,'body','hello','field2','world');

%%

%{
https://github.com/wjawaid/enrichR/blob/master/R/functions.R

        temp <- POST(url=paste0(getOption("enrichR.base.address"), "enrich"),
                     body=list(list=paste(genes, collapse="\n")))
%}
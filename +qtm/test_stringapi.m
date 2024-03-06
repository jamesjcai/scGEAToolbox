
% https://string-db.org/cgi/help?sessionId=bC9BtcIz7NfC
g = ["INS1","APOE","H2-K1","PLP1","CRYAB","S100A6","H2-D1","GATM","DBI","VIM","FXYD1","LY6A","RARRES2","HBB-BS","LDHB"];

s = urlencode(strtrim(sprintf('%s\r',g)));
% https://string-db.org/api/image/network?identifiers=DRD1_HUMAN%0dDRD2_HUMAN&species=9606
url=sprintf('https://string-db.org/api/image/network?identifiers=%s&species=9606',s);
web(url)
imshow(webread(url))


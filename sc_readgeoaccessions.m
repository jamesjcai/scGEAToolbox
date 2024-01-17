function [sce] = sc_readgeoaccessions(acc)

% if length(strsplit(acc,{',',';',' '}))>1
% end

url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s', acc);
a = webread(url);
b = strsplit(a, '\n');
b = string(b');

x = b((contains(b, 'Samples (')));
x = regexp(x, '\([0-9]+\)', 'match');
%x=strrep(x,'(','');
%x=strrep(x,')','');
x = str2double(x);
idx = find(contains(b, 'GSM'));
c = b(idx);
c = string(regexp(c, '>GSM[0-9]+', 'match'));
c = extractAfter(c, 1);

d = b(idx+1);
d = string(regexp(d, '>.*</', 'match'));
d = extractAfter(d, 1);
d = strrep(d, '/<', '');

assert(x == length(c));
assert(x == length(d));


[sce] = pkg.pipeline_multisamplesmerge(c);

% https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSM3308545
% https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE117770
% https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&id=303308547
% https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=303308547
% https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=200117770
% https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html


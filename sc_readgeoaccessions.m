function [sce]=sc_readgeoaccessions(acc)

% if length(strsplit(acc,{',',';',' '}))>1    
% end

url=sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s',acc);
a=webread(url);
b=strsplit(a,'\n');
b=string(b');

x=b((contains(b,'Samples (')));
x=regexp(x,'\([0-9]+\)','match');
%x=strrep(x,'(','');
%x=strrep(x,')','');
x=str2num(x);
idx=find(contains(b,'GSM'));
c=b(idx);
c=string(regexp(c,'>GSM[0-9]+','match'));
c=extractAfter(c,1);

d=b(idx+1);
d=string(regexp(d,'>.*</','match'));
d=extractAfter(d,1);
d=strrep(d,'/<','');

assert(x==length(c));
assert(x==length(d));
c

[sce]=pkg.pipeline_multisamplesmerge(c);


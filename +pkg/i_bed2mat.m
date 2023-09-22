function [chrid, startn, endn] = i_bed2mat(input)
%input=["1:12377-456";"20:234-56744";"X:234-4560"];
p = cell2mat(strfind(input, ":"));
q = cell2mat(strfind(input, "-"));
chrid = extractBefore(input, p);
startn = extractBetween(input, p+1, q-1);
endn = extractAfter(input, q);

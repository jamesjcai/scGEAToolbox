function [sce] = e_readgeotxt(~)

acc = 'GSM4318799';
url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s', acc);
a = webread(url);
b = strsplit(a, '\n');
c = string(b(contains(b, acc)))';
c = c(contains(c, 'ftp'));

if length(c) ~= 1
    disp(url)
    error('Unknown error.');
end

c1 = c(contains(c, 'txt'));
if isempty(c1), error('C1'); end
f1 = i_setupfile(c1);
if isempty(f1), error('f1'); end

[X, g] = sc_readtsvfile(f1);
sce = SingleCellExperiment(X, g);

end


function f = i_setupfile(c)
try
    tmpd = tempdir;
    [x] = regexp(c(1), '<a href="ftp://(.*)">(ftp', 'match');
    x = string(textscan(x, '<a href="ftp://%s'));
    x = append("https://", extractBefore(x, strlength(x)-5));
    if ~(ismcc || isdeployed)
        %#exclude urldecode
        x = urldecode(x);
    else
        x = pkg.urldecoding(x);
    end
    files = gunzip(x, tmpd);
    f = files{1};
catch
    f = [];
end
end

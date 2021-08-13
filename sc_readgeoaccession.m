function [sce]=sc_readgeoaccession(acc)

url=sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s',acc);
a=webread(url);
b=strsplit(a,'\n');
c=string(b(contains(b,acc)))';
c=c(contains(c,'ftp'));
if length(c)~=3 && length(c)~=1
    disp(url)
    error('Unknown error.');
end

switch length(c)
    case 3
        c1=c(contains(c,'mtx'));
        if isempty(c1), error('C1'); end
        f1=i_setupfile(c1);
        if isempty(f1), error('f1'); end
        c2=c(contains(c,'genes'));
        if isempty(c2), c2=c(contains(c,'features')); end
        if isempty(c2), error('C2'); end
        f2=i_setupfile(c2);
        if isempty(f2), error('f2'); end
        [X,g]=sc_readmtxfile(f1,f2);
    case 1
        c1=c(contains(c,'txt'));
        if isempty(c1)
            c1=c(contains(c,'csv'));
            if isempty(c1)
                c1=c(contains(c,'tsv'));
                if isempty(c1)
                    error('C1');
                end
            end
        end
            f1=i_setupfile(c1);
            if isempty(f1), error('f1'); end
            [X,g]=sc_readtsvfile(f1);
        
end
sce=SingleCellExperiment(X,g);
end


function f=i_setupfile(c)    
    try
        tmpd=tempdir;
        [x]=regexp(c(1),'<a href="ftp://(.*)">(ftp','match');
        x=string(textscan(x,'<a href="ftp://%s'));
        x=append("https://", extractBefore(x,strlength(x)-5));
        x=urldecode(x);
        fprintf('Downloading %s\n',x)
        files=gunzip(x,tmpd);
        f=files{1};
    catch
        f=[];
    end
end

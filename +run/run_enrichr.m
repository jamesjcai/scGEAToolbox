function run_enrichr(genelist,genenum)
% Run Enrichr
%
% see also: RUN_GORILLA

if nargin<1, genelist=[]; end
if nargin<2, genenum=1000; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/Enrichr');
addpath(pth);

infile=fullfile(pth,'input_template.html');
outfile=fullfile(pth,'input_page.html');

fid=fopen(infile,'r');
a=textscan(fid,'%s','delimiter','\n');
fclose(fid);
a=string(a{1});

fid=fopen(outfile,'w');
idx=find(a=='<textarea name=list rows=10 id=text-area cols=63></textarea>');
pause(1)
fprintf(fid,'%s\n',a(1:idx-1));

fprintf(fid,'<textarea name=list rows=10 id=text-area cols=63>');
n=min([length(genelist) genenum-1]);
if ~isempty(genelist)
    fprintf(fid,'%s\n',genelist(1));    
    for k=2:n
        fprintf(fid,'%s\n',genelist(k));
    end
end
fprintf(fid,'</textarea>\n');

fprintf(fid,'%s\n',a(idx+1:end));
fclose(fid);
pause(1)
web(outfile,'-browser');
pause(1)


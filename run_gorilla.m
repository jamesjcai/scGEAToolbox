function run_gorilla(genelist)
% Run GOrilla
% GOrilla is a tool for identifying and visualizing enriched GO terms in ranked lists of genes
% 
% http://cbl-gorilla.cs.technion.ac.il/
% Refs:
% Eran Eden*, Roy Navon*, Israel Steinfeld, Doron Lipson and Zohar Yakhini. "GOrilla: A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists", BMC Bioinformatics 2009, 10:48.
% Eran Eden, Doron Lipson, Sivan Yogev, Zohar Yakhini. "Discovering Motifs in Ranked Lists of DNA sequences", PLoS Computational Biology, 3(3):e39, 2007. 
%
% see also: RUN_ENRICHR
if nargin<1, genelist=[]; end

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/GOrilla');
addpath(pth);

infile=fullfile(pth,'input_template.html');
outfile=fullfile(pth,'input_page.html');

fid=fopen(infile,'r');
a=textscan(fid,'%s','delimiter','\n');
fclose(fid);
a=string(a{1});


fid=fopen(outfile,'w');
idx=find(a=="<TEXTAREA name=target_set rows=6 cols=63></TEXTAREA>");
pause(1)
fprintf(fid,'%s\n',a(1:idx-1));

fprintf(fid,'<TEXTAREA name=target_set rows=6 cols=63>');
if ~isempty(genelist)
    fprintf(fid,'%s\n',genelist(1));
    for k=2:length(genelist)
        fprintf(fid,'%s\n',genelist(k));
    end
end
fprintf(fid,'</TEXTAREA>\n');
fprintf(fid,'%s\n',a(idx+1:end));
fclose(fid);
pause(1)
web(outfile,'-browser');
pause(1)


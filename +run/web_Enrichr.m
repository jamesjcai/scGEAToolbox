function web_Enrichr(genelist, genenum, bkglist, wkdir)
% Run Enrichr
%
% see also: RUN.WEB_GORILLA, RUN.WEB_STRING

if nargin < 1, genelist = []; end
if nargin < 2, genenum = 100; end
if nargin < 3, bkglist = []; end
if nargin < 4, wkdir = ''; end

if isempty(wkdir), wkdir = tempdir; end

if ~isempty(bkglist)
    infile = fullfile(wkdir, 'input_template_bkg.html');
else
    infile = fullfile(wkdir, 'input_template.html');
end

[~, b]=fileparts(tempname);
% fx = sprintf('input_page_%s.html', char(randi([97, 122], 1, 8)));
fx = sprintf('input_page_%s.html', b);
outfile = fullfile(wkdir, fx);

fid = fopen(infile, 'r');
a = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
a = string(a{1});

fid = fopen(outfile, 'w');
idx = find(a == '<textarea name=list rows=10 id=text-area cols=63></textarea>');
pause(1)
fprintf(fid, '%s\n', a(1:idx-1));

fprintf(fid, '<textarea name=list rows=10 id=text-area cols=63>');

n = min([length(genelist), genenum]);

% n=length(genelist);

if ~isempty(genelist)
    if isstring(genelist)
        fprintf(fid, '%s\n', genelist(1));
        for k = 2:n
            fprintf(fid, '%s\n', genelist(k));
        end
    elseif iscell(genelist)
        fprintf(fid, '%s\n', genelist{1});
        for k = 2:n
            fprintf(fid, '%s\n', genelist{k});
        end
    end
end
fprintf(fid, '</textarea>\n');

fprintf(fid, '%s\n', a(idx+1:end));
fclose(fid);
pause(1)
web(outfile, '-browser');
pause(1)

end
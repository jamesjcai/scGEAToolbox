function web_Enrichr_bkg(genelist, bkglist, genenum, wkdir)
% Run Enrichr
%
% see also: RUN.WEB_GORILLA, RUN.WEB_STRING

if nargin < 1, genelist = []; end
if nargin < 2, bkglist = []; end
if nargin < 3, genenum = 200; end
if nargin < 4, wkdir = ''; end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'external', 'web_Enrichr');

if ~isempty(bkglist)
    infile = fullfile(pth, 'input_template_bkg.html');
else
    infile = fullfile(pth, 'input_template.html');
end


[~, b]=fileparts(tempname);
% fx = sprintf('input_page_%s.html', char(randi([97, 122], 1, 8)));
fx = sprintf('input_page_%s.html', b);
if ~isempty(wkdir)
    outfile = fullfile(wkdir, fx);
else    
    outfile = fullfile(pth, fx);
end


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

b=a(idx+1:end);


idx = find(b == '<textarea name=background rows=15 id=text-area cols=63></textarea>');
pause(1)
fprintf(fid, '%s\n', b(1:idx-1));

fprintf(fid, '<textarea name=background rows=15 id=text-area cols=63>');

% n = min([length(genelist), genenum]);

n=length(bkglist);

if ~isempty(bkglist)
    if isstring(bkglist)
        fprintf(fid, '%s\n', bkglist(1));
        for k = 2:n
            fprintf(fid, '%s\n', bkglist(k));
        end
    elseif iscell(bkglist)
        fprintf(fid, '%s\n', bkglist{1});
        for k = 2:n
            fprintf(fid, '%s\n', bkglist{k});
        end
    end
end
fprintf(fid, '</textarea>\n');

fprintf(fid, '%s\n', b(idx+1:end));
fclose(fid);
pause(1)
web(outfile, '-browser');
pause(1)

end
function web_STRING(genelist,genenum,species)
% Run STRING
%
% see also: RUN.WEB_GORILLA, RUN.WEB_ENRICHR


pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'external','web_STRING');
%infile=fullfile(pth,'input_template.html');

if nargin<1, genelist=[]; end
if nargin<2, genenum=50; end
if nargin<3, species="human"; end

switch lower(species)
    case 'human'
        species="9606";
    case 'mouse'
        species="10090";
    case 'rat'
        species="10116";
    case 'zebrafish'
        species="10116";
    otherwise
        taxonid = pkg.i_species2taxid(species);
        species=sprintf("%d",taxonid);
end

n=min([length(genelist) genenum-1]);
genelist=genelist(1:n);
if ~isempty(genelist)
    s=sprintf('%s%%0d',genelist);
end

url=sprintf('https://string-db.org/api/svg/network?identifiers=%s&species=%s',s,species);
web(url);
% ref: https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
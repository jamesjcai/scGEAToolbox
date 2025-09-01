function [r] = sc_potency(X, genelist, speciesid)
%Estimate differentiation potency of cells
%without the need to assume prior biological knowledge
%such as marker expression or timepoint.
% CCAT (Correlation of Connectome And Transcriptome)
%https://github.com/aet21/SCENT
%
%see also: RUN.CYTOTRACE

if nargin < 3, speciesid = 1; end

if ~isnumeric(speciesid)
    switch speciesid
        case 'human'
            speciesid = 1;
        case 'mouse'
            speciesid = 2;
        otherwise
            speciesid = 1;
    end
end

if ~ismember(speciesid, [1, 2])
    error('Invalid SPECIESID');
end

pw1 = fileparts(mfilename('fullpath'));


dbfile1 = fullfile(pw1, 'assets', 'STRING', 'stringdb_human.mat');
dbfile2 = fullfile(pw1, 'assets', 'STRING', 'stringdb_mouse.mat');

if ~exist(dbfile1, 'file')
    if ~exist(fileparts(dbfile1), 'dir')
        mkdir(fileparts(dbfile1));
    end
    %disp('Downloading ...... stringdb_human.mat')
    url = 'https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/stringdb_human.mat';
    % outfilename =
    websave(dbfile1, url);
end
if ~exist(dbfile2, 'file')
    if ~exist(fileparts(dbfile2), 'dir')
        mkdir(fileparts(dbfile2));
    end
    %disp('Downloading ...... stringdb_mouse.mat')
    url = 'https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/stringdb_mouse.mat';
    %outfilename =
    websave(dbfile2, url);
end

genelist = upper(genelist);
if speciesid == 1
    ppinetfile = dbfile1; % 'Z:\Cailab\CCC_utilities\STRING\stringdb_human.mat';
    disp('Using stringdb_human.mat')
else
    ppinetfile = dbfile2; % 'Z:\Cailab\CCC_utilities\STRING\stringdb_mouse.mat';
    disp('Using stringdb_mouse.mat')
end
load(ppinetfile, 'G');
G.Edges.Weight = double(G.Edges.Weight > 0);
GNodes = upper(string(table2array(G.Nodes)));
Gdegree = G.degree;

if issparse(X), X = full(X); end
X = log2(X+1.1);

[~, i, j] = intersect(genelist, GNodes);
if isempty(i)
    error(sprintf('No gene is found in STRING PPI data base.\nCheck gene names.'));
end
d = Gdegree(j);
X = X(i, :);
r = corr(X, d); % Correlation of Connectome And Transcriptome
end

function [T] = e_genecoordinates(genome)
%E_GENECOORDINATES Genomic coordinates of genes for a reference genome.
%
%   T = PKG.E_GENECOORDINATES(genome) returns a table with one row per gene
%   and the variables:
%
%       Gene        gene symbol (falls back to the Ensembl ID when the
%                   Ensembl release carries no symbol for that gene)
%       EnsemblID   Ensembl gene stable ID
%       Chr         chromosome number (1..22 human, 1..19 mouse)
%       Start       start coordinate, bp
%       Stop        stop coordinate, bp
%
%   Rows are sorted by Chr and then by Start, which is the order a CNV
%   moving average along the genome expects. Only autosomes are covered.
%
%   genome is 'hg38' (default), 'hg19' or 'mm10'.
%
%   The table is read from assets/scCancer/gene_chr_<genome>.txt and cached,
%   so repeated calls within a session do not re-read the file.
%
%   See also SC_INFERCNV, PKG.E_ENSEMBL2SYMBOL.

if nargin < 1 || isempty(genome), genome = 'hg38'; end
genome = lower(string(genome));

persistent cachedGenome cachedTable
if ~isempty(cachedGenome) && cachedGenome == genome
    T = cachedTable;
    return;
end

switch genome
    case {"hg38", "hg19"}
        species = 'human';
    case "mm10"
        species = 'mouse';
    otherwise
        error("Unknown genome %s. Use 'hg38', 'hg19' or 'mm10'.", genome);
end

pw1 = fileparts(mfilename('fullpath'));
dbfile = fullfile(pw1, '..', 'assets', 'scCancer', ...
    sprintf('gene_chr_%s.txt', genome));
if ~exist(dbfile, 'file')
    error('Missing file %s.', dbfile);
end

T = readtable(dbfile, 'FileType', 'text', 'Delimiter', '\t', ...
    'ReadVariableNames', false, 'TextType', 'string');
T.Properties.VariableNames = {'EnsemblID', 'Chr', 'Start', 'Stop'};

T.Chr = double(extractAfter(T.Chr, 3));
T.Gene = pkg.e_ensembl2symbol(T.EnsemblID, species);

% e_ensembl2symbol leaves IDs it cannot map untouched, so a row whose Gene
% still equals its EnsemblID simply has no symbol in this Ensembl release.
T = T(:, {'Gene', 'EnsemblID', 'Chr', 'Start', 'Stop'});
T = sortrows(T, {'Chr', 'Start'});

cachedGenome = genome;
cachedTable = T;
end

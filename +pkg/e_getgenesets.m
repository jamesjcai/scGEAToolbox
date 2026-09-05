function [setmatrx, setnames, setgenes] = e_getgenesets(option,species,parentfig,confidence)
% E_GETGENESETS  Load a gene set collection as a membership matrix.
%
%   [setmatrx, setnames, setgenes] = E_GETGENESETS(option, species, parentfig)
%   returns an nSets-by-nGenes membership matrix, the set names, and the gene
%   symbols labelling the columns.
%
%   ... = E_GETGENESETS(option, species, parentfig, confidence) restricts the
%   DoRothEA options to the given confidence levels, e.g. "ABC", 'A', or
%   ["A" "B"]. Omit or pass [] to keep every level, which is the default and
%   the historical behaviour. DoRothEA is heavily weighted toward its lowest
%   level -- A 5434 edges, B 974, C 6815, D 15863, E 425418 -- so the whole
%   collection is 94% level E. "ABC" is the decoupleR convention for regulon
%   activity inference. CONFIDENCE is ignored, with a warning, for the
%   non-DoRothEA options.
%
%   OPTION:
%     1 / 'MSIGDB'    - MSigDB molecular signatures (opens a selector dialog)
%     2 / 'TF'        - DoRothEA TF targets, unsigned: activated edges only,
%                       returned logical
%     3 / 'Predefined'- scGEAToolbox predefined marker programs
%     4 / 'Glycobiology' - curated glycobiology gene sets
%     5 / 'TFsigned'  - DoRothEA TF targets, signed: +1 activated and -1
%                       repressed, returned sparse double. Use this with
%                       methods that accept weighted set membership
%                       (SC_GSETTEST Method="ulm", PKG.E_ULM); the unsigned
%                       option 2 cannot detect a regulon whose targets move
%                       in opposite directions.
%
% See also: SC_GSETTEST, SC_DPG, PKG.E_ULM, PKG.E_GLYCOGENESETS
if nargin < 4, confidence = []; end
if nargin < 3, parentfig = []; end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end
if nargin<2, species=[]; end
if isempty(species) && ~isequal(option, 1) && ~isequal(option, 'MSIGDB') && ~isequal(option, 'MSigDB Molecular Signatures')
    species = 'human';
end
if nargin<1 || isempty(option), option = 1; end

conflevels = i_normalizeconfidence(confidence);
if ~isempty(conflevels) && ~i_isdorotheaoption(option)
    warning("pkg:e_getgenesets:ConfidenceIgnored", ...
        "CONFIDENCE applies only to the DoRothEA options " + ...
        "(2/'TF' and 5/'TFsigned'); it is ignored for option '%s'.", ...
        string(option));
end

setmatrx=[];
setnames=[];
setgenes=[];
switch option
    case {1,'MSIGDB','MSigDB Molecular Signatures'}
        % [Col]=pkg.e_getmsigdbset;
        [~, ~, Col, ctag] = gui.i_selectMSigDBGeneSet(species,true,parentfig);
        if isempty(Col) || isempty(ctag)
            return;
        end
        pw1 = fileparts(mfilename('fullpath'));
        isloaded=false;
        try
            dbfile = fullfile(pw1, '..', 'assets', 'MSigDB', ...
                sprintf('msigdb_%s.mat',ctag));
            load(dbfile,'setmatrx','setnames','setgenes');
            isloaded=true;
        catch ME
            warning(ME.message);
        end

        if ~isloaded
            setnames = string(fields(Col));
            gcell = cell(length(setnames), 1);
            for k=1:length(setnames)
                gcell{k} = Col.(setnames{k}).geneSymbols;
            end
            glist = unique(vertcat(gcell{:}));
            setgenes=glist(strlength(glist)>0);
            setmatrx=false(length(setnames),length(setgenes));
            for k=1:length(setnames)
                tgsPos = string(Col.(setnames(k)).geneSymbols);
                setmatrx(k,:)=ismember(setgenes,tgsPos);
            end
            % save(sprintf('msigdb_%s',ctag),'setmatrx','setnames','setgenes');
        end

    case {'TF',2,'DoRothEA TF Targets'}
        T = i_loaddorothea(conflevels);
        Ttfgn = T(T.mor > 0, :);
        [gid, gnlist] = findgroups(string(Ttfgn.target));
        [tid, tflist] = findgroups(string(Ttfgn.tf));
        t = zeros(max(tid), max(gid));
        t(sub2ind([max(tid), max(gid)], tid, gid)) = Ttfgn.mor;

        setmatrx=logical(t);
        setnames=string(tflist);
        setgenes=string(gnlist);
    case {5,'TFsigned','DoRothEA TF Targets (signed)'}
        % Signed regulons. Unlike option 2 this keeps the repressing edges
        % and preserves the mode of regulation. Returned sparse because the
        % dense form is 1333-by-20295 doubles, about 216 MB.
        T = i_loaddorothea(conflevels);
        [gid, gnlist] = findgroups(string(T.target));
        [tid, tflist] = findgroups(string(T.tf));
        setmatrx = sparse(tid, gid, double(T.mor), numel(tflist), numel(gnlist));
        setnames = string(tflist);
        setgenes = string(gnlist);
    case {3,'Predefined'}
        [~, T] = pkg.e_cellscores([], [], 0);
        gcell2 = cell(height(T), 1);
        for k=1:height(T)
            gcell2{k} = unique(strsplit(string(T.PositiveMarkers(k)), ','))';
        end
        glist = unique(vertcat(gcell2{:}));
        setgenes=glist(strlength(glist)>0);
        setnames=string(T.ScoreType);
        setmatrx=false(length(setnames),length(setgenes));
        for k=1:length(setnames)
            tgsPos = unique(strsplit(string(T.PositiveMarkers(k)), ','));
            setmatrx(k,:)=ismember(setgenes,tgsPos);
        end
    case {4,'Glycobiology','Glycobiology Gene Sets'}
        pw1 = fileparts(mfilename('fullpath'));
        glycofile = fullfile(pw1, '..', 'assets', 'GeneSets', ...
            'glycobiology.mat');
        if isfile(glycofile)
            load(glycofile, 'setmatrx', 'setnames', 'setgenes');
        else
            % Fall back to building the collection on the fly.
            [setmatrx, setnames, setgenes] = pkg.e_glycogenesets();
        end
end
end

% =========================================================================
function T = i_loaddorothea(conflevels)
% Load the DoRothEA TF-target table, optionally restricted by confidence.
pw1 = fileparts(mfilename('fullpath'));
fname = fullfile(pw1, '..', 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
S = load(fname, 'T');
T = S.T;
if isempty(conflevels)
    return;
end
have = unique(string(T.confidence));
unknown = setdiff(conflevels, have);
if ~isempty(unknown)
    error("pkg:e_getgenesets:BadConfidence", ...
        "Unknown DoRothEA confidence level(s) %s. Available levels are %s.", ...
        join(unknown, ", "), join(sort(have), ", "));
end
T = T(ismember(string(T.confidence), conflevels), :);
end

% =========================================================================
function levels = i_normalizeconfidence(confidence)
% Accept "ABC", 'A', ["A" "B"] or [] and return unique single-letter levels.
if isempty(confidence)
    levels = strings(0, 1);
    return;
end
c = erase(upper(string(confidence(:))), [" ", ","]);
levels = unique(string(num2cell(char(join(c, "")))));
levels = levels(strlength(levels) > 0);
levels = levels(:);
end

% =========================================================================
function tf = i_isdorotheaoption(option)
tf = isequal(option, 2) || isequal(option, 5) || ...
    any(strcmpi(string(option), ["TF", "TFsigned", ...
    "DoRothEA TF Targets", "DoRothEA TF Targets (signed)"]));
end

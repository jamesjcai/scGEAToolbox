function [setmatrx, setnames, setgenes] = e_getgenesets(option,species)
if nargin<2 || isempty(species), species='human'; end
if nargin<1 || isempty(option), option = 1; end

setmatrx=[];
setnames=[];
setgenes=[];
switch option
    case {1,'MSIGDB','MSigDB Molecular Signatures'}
        % [Col]=pkg.e_getmsigdbset;
        [~, ~, Col, ctag] = gui.i_selectMSigDBGeneSet(species,true);
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
        pw1 = fileparts(mfilename('fullpath'));
        fname = fullfile(pw1, '..','assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
        load(fname, 'T');
        Ttfgn = T(T.mor > 0, :);
        [gid, gnlist] = findgroups(string(Ttfgn.target));
        [tid, tflist] = findgroups(string(Ttfgn.tf));
        t = zeros(max(tid), max(gid));
        t(sub2ind([max(tid), max(gid)], tid, gid)) = Ttfgn.mor;

        setmatrx=logical(t);
        setnames=string(tflist);
        setgenes=string(gnlist);
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
end

function [setmatrx, setnames, setgenes] = e_getgenesets(option,species)
if nargin<2 || isempty(species), species='human'; end
if nargin<1 || isempty(option), option = 1; end

setmatrx=[];
setnames=[];
setgenes=[];
switch option
    case {1,'MSIGDB'}
        % [Col]=pkg.e_getmsigdbset;
        [~, ~, Col, ctag] = gui.i_selectMSigDBGeneSet(species,true);
        if isempty(Col) || isempty(ctag)
            return;
        end

        pw1 = fileparts(mfilename('fullpath'));
        isloaded=false;
        try
            dbfile = fullfile(pw1, '..', 'resources', 'MSigDB', ...
                sprintf('msigdb_%s.mat',ctag));
            load(dbfile,'setmatrx','setnames','setgenes');
            isloaded=true;
        catch ME
            warning(ME.message);            
        end

        if ~isloaded
            setnames = string(fields(Col));
            glist=[];
            for k=1:length(setnames)  
                glist=[glist;Col.(setnames{k}).geneSymbols];
            end
            glist=unique(glist);
            setgenes=glist(strlength(glist)>0);
            setmatrx=false(length(setnames),length(setgenes));
            for k=1:length(setnames)
                tgsPos = string(Col.(setnames(k)).geneSymbols);
                setmatrx(k,:)=ismember(setgenes,tgsPos);
            end
            % save(sprintf('msigdb_%s',ctag),'setmatrx','setnames','setgenes');
        end
   
    case {'TF',2}
        pw1 = fileparts(mfilename('fullpath'));
        fname = fullfile(pw1, '..','resources', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
        load(fname, 'T');
        Ttfgn = T(T.mor > 0, :);
        [gid, gnlist] = grp2idx(Ttfgn.target);
        [tid, tflist] = grp2idx(Ttfgn.tf);
        t = zeros(max(tid), max(gid));
        t(sub2ind([max(tid), max(gid)], tid, gid)) = Ttfgn.mor;
        
        setmatrx=logical(t);
        setnames=string(tflist);
        setgenes=string(gnlist);
    case {3,'Predefined'}
        [~, T] = pkg.e_cellscores([], [], 0);
        glist=[];
        for k=1:size(T,1)
            tgsPos = unique(strsplit(string(T.PositiveMarkers(k)), ','));
            glist=[glist;tgsPos'];
        end
        glist=unique(glist);
        setgenes=glist(strlength(glist)>0);
        setnames=string(T.ScoreType);
        setmatrx=false(length(setnames),length(setgenes));
        for k=1:length(setnames)
            tgsPos = unique(strsplit(string(T.PositiveMarkers(k)), ','));
            setmatrx(k,:)=ismember(setgenes,tgsPos);
        end
end



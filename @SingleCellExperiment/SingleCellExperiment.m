classdef SingleCellExperiment
   properties
      X double % {mustBeNumeric, mustBeFinite, mustBeNonNan} % counts
      g string                               % genelist
      s double {mustBeNumeric, mustBeFinite} % cell embeddings
      c                                      % current/active group/class id
      c_cell_cycle_tx                        % cell cycle string
      c_cell_type_tx                         % cell type string
      c_cluster_id                           % clustering result
      c_batch_id                             % batch id
      c_cell_id                              % barcode
      list_cell_attributes cell  % e.g., attributes = {'size',[4,6,2]};
      list_gene_attributes cell  % e.g., attributes = {'size',[4,6,2]};
      metadata string
      struct_cell_embeddings=struct('tsne',[],'umap',[], ...
          'phate',[],'metaviz',[])
      struct_cell_clusterings=struct('kmeans',[],'snndpc',[],...
                                      'sc3',[],'simlr',[],'soptsc',[],...
                                      'sinnlrr',[],'specter',[],...
                                      'seurat',[])
      table_attributes table      
   end
   
   % inputdlg('Sample Info:','Meta data',[10 50],{char(b.')})

   properties (Dependent)
      NumCells
      NumGenes
   end 

   methods

   function obj = SingleCellExperiment(X,g,s,c)
        if nargin<1, X=[]; end
        if nargin<2 || isempty(g), g=pkg.i_num2strcell(size(X,1),"g"); end
        if nargin<3 || isempty(s), s=randn(size(X,2),3); end
        if nargin<4 || isempty(c), c=ones(size(X,2),1); end
        assert(size(X,2)==size(s,1))
        if ~(size(s,2)>2)
            s=[s zeros(size(X,2),1)];
        end
        obj.X = X;
        obj.g=g;
        obj.s=s;
        obj.c=c;
        obj.c_cell_id=transpose(1:size(X,2));
        obj.c_batch_id=ones(size(X,2),1);
        obj.c_cluster_id=ones(size(X,2),1);
        obj.c_cell_cycle_tx=repmat("undetermined",size(X,2),1);
        obj.c_cell_type_tx=repmat("undetermined",size(X,2),1);
        obj.metadata=string(sprintf('Created: %s',datetime()));
        if ~isfield(obj.struct_cell_embeddings,'metaviz')
            obj.struct_cell_embeddings.('metaviz')=[];
        end
        % obj.struct_cell_embeddings=struct('tsne',[],'umap',[],'phate',[]);
    end

   function m = get.NumCells(obj)
      m = size(obj.X,2); 
   end
   
   function obj = set.NumCells(obj,~)
      fprintf('%s%d\n','NumCells is: ',obj.NumCells)
      error('You cannot set NumCells property'); 
   end
   
   function m = get.NumGenes(obj)   
      m = size(obj.X,1); 
   end
   
   function obj = set.NumGenes(obj,~)
      fprintf('%s%d\n','NumGenes is: ',obj.NumGenes)
      error('You cannot set NumGenes property'); 
   end
 
    function r=numcells(obj)
        r=size(obj.X,2);
    end
    
    function r=numgenes(obj)
        r=size(obj.X,1);
    end
    
    obj = estimatepotency(obj,speciesid,forced)
    obj = estimatecellcycle(obj,forced,methodid)
    obj = embedcells(obj,methodid,forced,usehvgs,ndim,numhvg,whitelist)
    obj = clustercells(obj,k,methodid,forced)
    obj = assigncelltype(obj,speciesid)
    obj = qcfilterwhitelist(obj,libszcutoff,mtratio,min_cells_nonzero,gnnumcutoff,whitelist)
    obj = sortcells(obj,idx);
    
    function obj = removecells(obj,idx)
            obj.X(:,idx)=[];
            obj.s(idx,:)=[];
            obj.c(idx)=[];
            if ~isempty(obj.c_cell_cycle_tx)
                obj.c_cell_cycle_tx(idx)=[];
            end
            if ~isempty(obj.c_cell_type_tx)
                obj.c_cell_type_tx(idx)=[];
            end
            if ~isempty(obj.c_cluster_id)
                obj.c_cluster_id(idx)=[];
            end
            if ~isempty(obj.c_batch_id)
                obj.c_batch_id(idx)=[];
            end
            if ~isempty(obj.c_cell_id)
                obj.c_cell_id(idx)=[];
            end
            for k=2:2:length(obj.list_cell_attributes)
                obj.list_cell_attributes{k}(idx)=[];
            end
            
            a=fields(obj.struct_cell_embeddings);
            for k=1:length(a)
                if ~isempty(obj.struct_cell_embeddings.(a{k}))
                    obj.struct_cell_embeddings.(a{k})(idx,:)=[];
                end
            end
            
            a=fields(obj.struct_cell_clusterings);
            for k=1:length(a)
                if ~isempty(obj.struct_cell_clusterings.(a{k}))
                     obj.struct_cell_clusterings.(a{k})(idx)=[];
                end
            end
            % obj.NumCells=size(obj.X,2);
            if ~isempty(obj.table_attributes)
                if istable(obj.table_attributes)
                    if size(obj.table_attributes,1)==obj.NumCells
                        obj.table_attributes(idx,:)=[];
                    end
                end
            end
        end

    function obj = selectcells(obj,idx)
        if islogical(idx) && length(idx)==obj.NumCells
            ix=idx;
        else
            ix=false(obj.NumCells,1);
            ix(idx)=true;            
        end
        obj = removecells(obj,~ix);
    end

    function obj = set.c(obj,tmpc)
      if length(tmpc)~=numcells(obj)
         error('length(c)~=numcells(sce)');
      else
         obj.c=tmpc;
      end
    end

    function r = title(obj)
%        r=sprintf('%d x %d\n[genes x cells]',...
%            size(obj.X,1),size(obj.X,2));

        r=sprintf('%d x %d',...
           size(obj.X,1),size(obj.X,2));
    end
    
    function obj = qcfilter(obj,libsize,mtratio,min_cells_nonzero)
        if nargin<4 || isempty(min_cells_nonzero), min_cells_nonzero=15; end
        if nargin<3 || isempty(mtratio), mtratio=0.15; end
        if nargin<2 || isempty(libsize), libsize=1000; end
        %        case 'Relaxed (keep more cells/genes)'
        %            definput = {'500','0.20','10'};
        %        case 'Strigent (keep less cells/genes)'
        %            definput = {'1000','0.15','15'};
        [~,keptg,keptidxv]=sc_qcfilter(obj.X,obj.g,libsize,mtratio,...
                                       min_cells_nonzero);
        for k=1:length(keptidxv)
            obj = selectcells(obj,keptidxv{k});
        end
        [y]=ismember(obj.g,keptg);
        obj.X=obj.X(y,:);
        obj.g=obj.g(y);
        [obj.X,obj.g]=sc_rmdugenes(obj.X,obj.g);
    end
    
    function obj = selectgenes(obj,min_cellnum,nonzero_cutoff)
        if nargin<3, nonzero_cutoff=1; end
        if nargin<2, min_cellnum=0.01; end
        [tmpX,tmpg]=sc_selectg(obj.X,obj.g,min_cellnum,nonzero_cutoff);
        obj.X=tmpX;
        obj.g=tmpg;
    end
    
    function obj = selectkeepgenes(obj,min_countnum,min_cellnum)
        if nargin<2, min_countnum=1; end
        if nargin<3, min_cellnum=0.01; end        
        nc=sum(obj.X>=min_countnum,2);
        if min_cellnum<1
            idxkeep1=nc>=min_cellnum*size(obj.X,2);
        else
            idxkeep1=nc>=min_cellnum;
        end
        idxkeep2=true(size(idxkeep1));
        k=sum(~idxkeep1);
        [~,idx2]=mink(mean(obj.X,2),k);
        idxkeep2(idx2)=false;
        idxkeep=idxkeep1|idxkeep2;
        obj.X=obj.X(idxkeep,:);
        obj.g=obj.g(idxkeep);
    end
    
    function obj = rmmtgenes(obj)
        [tmpX,tmpg,idx]=sc_rmmtgenes(obj.X,obj.g,'mt-',true);
        if sum(idx)>0
            obj.X=tmpX;
            obj.g=tmpg;
        end
    end
    
    function obj = rmribosomalgenes(obj)
        ribog=pkg.i_get_ribosomalgenes;
        [i]=ismember(upper(obj.g),ribog);        
        obj.X=obj.X(~i,:);
        obj.g=obj.g(~i);
        fprintf('%d ribosomal genes found and removed.\n',...
            sum(i));
    end

    function obj = appendmetainfo(obj,infostr)
        if ~isstring(infostr)
            infostr=string(infostr);
        end
        obj.metadata = [obj.metadata; infostr];
    end

    function c_check( self )
        assert(~isempty( self.c ),'Must be defined!');
    end
    
    % function disp(td)
    %   fprintf(1,...
    %      'SingleCellExperiment: %d genes x %d cells\n',...
    %      numgenes(td),numcells(td));
    % end 
    
   end  

% https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html   
end



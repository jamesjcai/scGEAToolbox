classdef SingleCellExperiment
   properties
      X double {mustBeNumeric, mustBeFinite} % counts
      g string                               % genelist
      s double {mustBeNumeric, mustBeFinite} % cell.embeddings
      c                                      % current group/class id
      c_cell_cycle_tx                        % cell cycle string
      c_cell_type_tx                         % cell type string
      c_cluster_id                           % clustering result
      c_batch_id                             % batch id
      c_cell_id                              % barcode
      list_cell_attributes cell  % e.g., attributes = {'size',[4,6,2]};
      list_gene_attributes cell  % e.g., attributes = {'size',[4,6,2]};
      struct_cell_embeddings=struct('tsne',[],'umap',[],'phate',[])
      struct_cell_clusterings=struct('kmeans',[],'snndpc',[],...
                                      'sc3',[],'simlr',[],'soptsc',[],...
                                      'sinnlrr',[],'specter',[])
      table_attributes table
   end
   
   properties (Dependent)
      NumCells
      NumGenes
   end 

   methods
   % output = myFunc(obj,arg1,arg2)
   function obj = SingleCellExperiment(X,g,s,c)
        if nargin<1, X=[]; end
        if nargin<2 || isempty(g), g=string(transpose(1:size(X,1))); end
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
    obj = embedcells(obj,methodid,forced,usehvgs)
    obj = clustercells(obj,k,methodid,forced)
    obj = assigncelltype(obj,speciesid)
    
    function obj = removecells(obj,i)
            obj.X(:,i)=[];
            obj.s(i,:)=[];
            obj.c(i)=[];
            if ~isempty(obj.c_cell_cycle_tx)
                obj.c_cell_cycle_tx(i)=[];
            end
            if ~isempty(obj.c_cell_type_tx)
                obj.c_cell_type_tx(i)=[];
            end
            if ~isempty(obj.c_cluster_id)
                obj.c_cluster_id(i)=[];
            end
            if ~isempty(obj.c_batch_id)
                obj.c_batch_id(i)=[];
            end
            if ~isempty(obj.c_cell_id)
                obj.c_cell_id(i)=[];
            end
            for k=2:2:length(obj.list_cell_attributes)
                obj.list_cell_attributes{k}(i)=[];
            end
        end    

    function obj = selectcells(obj,i)
        if islogical(i) && length(i)==obj.NumCells
            ix=i;
        else
            ix=true(obj.NumCells,1);
            ix(i)=false;            
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
       r=sprintf('%d x %d\n[genes x cells]',...
           size(obj.X,1),size(obj.X,2));
    end
    
    function obj = qcfilter(obj)
        [~,keptg,keptidxv]=sc_qcfilter(obj.X,obj.g);
        for k=1:length(keptidxv)
            obj = selectcells(obj,keptidxv{k});
        end
        [y]=ismember(obj.g,keptg);
        obj.X=obj.X(y,:);
        obj.g=obj.g(y);
        [obj.X,obj.g]=sc_rmdugenes(obj.X,obj.g);
    end
    
    function obj = selectgenes(obj,min_countnum,min_cellnum)
        if nargin<2, min_countnum=1; end
        if nargin<3, min_cellnum=0.01; end
        [tmpX,tmpg]=sc_selectg(obj.X,obj.g,min_countnum,min_cellnum);
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



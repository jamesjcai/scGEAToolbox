classdef SingleCellExperiment
   properties
      X double {mustBeNumeric, mustBeFinite} % counts
      g string                               % genelist
      s double {mustBeNumeric, mustBeFinite} % cell.embeddings
      c                                      % cell group/class id
      c_cell_cycle_phase_tx
      c_cell_type_tx
      c_cluster_id
      c_batch_id
      c_cell_id
      list_cell_attributes cell  % e.g., attributes = {'size',[4,6,2]};
      list_gene_attributes cell  % e.g., attributes = {'size',[4,6,2]};
      table_attributes table
      % drmethod {mustBeMember(drmethod,{'tsne','umap','phate'})} = 'tsne'
   end
   
   properties (Dependent)
      NumCells
      NumGenes
   end 

   methods
    function obj = SingleCellExperiment(X,g,s,c)
        if nargin<1, X=[]; end
        if nargin<2 || isempty(g), g=string(transpose(1:size(X,1))); end
        if nargin<3 || isempty(s), s=randn(size(X,2),3); end
        if nargin<4 || isempty(c), c=ones(size(X,2),1); end
        obj.X = X;
        obj.g=g;
        obj.s=s;
        obj.c=c;
        obj.c_cell_id=transpose(1:size(X,2));
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
    
    function obj = estimatepotency(obj)
        if sum(strcmp('cell_potency',obj.list_cell_attributes))==0;
            idx=input('Species: 1=human,2=mouse >>');
            r=sc_potency(obj.X,obj.g,idx);
            obj.list_cell_attributes=[obj.list_cell_attributes,...
                {'cell_potency',r}];
            disp('cell_potency added.');
        else
            disp('cell_potency existed.');
        end
    end

    function obj = estimatecellcycle(obj,forced)
        if nargin<2, forced=false; end
        if isempty(obj.c_cell_cycle_phase_tx) || forced
            obj.c_cell_cycle_phase_tx=run_cellcycle(obj.X,obj.g);
            disp('c_cell_cycle_phase_tx added.');
        else
            disp('c_cell_cycle_phase_tx existed.');
        end
    end
    
    function obj = removecells(obj,i)
            obj.X(:,i)=[];
            obj.s(i,:)=[];
            obj.c(i)=[];    
            if ~isempty(obj.c_cell_cycle_phase_tx)
                obj.c_cell_cycle_phase_tx(i)=[];
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

    function r=title(obj)
       r=sprintf('%d x %d\n[genes x cells]',...
           size(obj.X,1),size(obj.X,2));
    end
    
    function obj = qcfilter(obj)
        [~,~,keptidxv]=sc_qcfilter(obj.X,obj.g);
        for k=1:length(keptidxv)
            obj = selectcells(obj,keptidxv{k});
        end
    end
    
    function obj = selectgenes(obj,min_countnum,min_cellnum)
        if nargin<2, min_countnum=1; end
        if nargin<3, min_cellnum=0.05; end
        [tmpX,tmpg]=sc_selectg(obj.X,obj.g,min_countnum,min_cellnum);
        obj.X=tmpX;
        obj.g=tmpg;
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

    
% function disp(td)
%   fprintf(1,...
%      'SingleCellExperiment: %d genes x %d cells\n',...
%      numgenes(td),numcells(td));
% end 
   end

% https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html   
end



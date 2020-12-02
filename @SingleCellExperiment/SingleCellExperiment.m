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

    function r = libsz(obj)
     r = sum([obj.X]);
    end
    function r = norm(obj)
     r = sc_norm(obj.X);
    end
    
    function r=numcells(obj)
        r=size(obj.X,2);
    end    
    function r=numgenes(obj)
        r=size(obj.X,1);
    end

    function obj = rmcells(obj,i)
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
    end

    function obj = selectcells(obj,i)
        obj = rmcells(obj,~i);
    end
    
   function obj = set.c(obj,cx)
      if length(cx)~=numcells(obj)         
         error('You cannot set the Modulus property');
      else
         obj.c=cx;
      end
   end 
   function r=title(obj)
       r=sprintf('%d x %d\n[genes x cells]',...
           size(obj.X,1),size(obj.X,2));
   end
   
% function disp(td)
%   fprintf(1,...
%      'SingleCellExperiment: %d genes x %d cells\n',...
%      numgenes(td),numcells(td));
% end
   
end
% https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html   
end


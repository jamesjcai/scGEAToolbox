classdef SingleCellExperiment
   properties
      X double {mustBeNumeric, mustBeFinite}
      genelist 
      drmethod {mustBeMember(drmethod,{'tsne','umap','phate'})} = 'tsne'
      s (:,3)
      c
      
      c_cell_cycle_phase
      c_cell_type
      c_cluster_id
   end
      properties (Dependent)
      c_cell_idx
   end 
   methods
    function obj = SingleCellExperiment(X,genelist,s,c)
        if nargin == 1
            obj.X = X;
            obj.genelist=string((1:size(X,1))');
            obj.s=randn(size(X,2),3);
            obj.c=ones(size(X,2),1);
        elseif nargin==2
            obj.X = X;
            obj.genelist=genelist;
            obj.s=randn(size(X,2),3);
            obj.c=ones(size(X,2),1);
        elseif nargin==3
            obj.X = X;
            obj.genelist=genelist;
            obj.s=s;
            obj.c=ones(size(X,2),1);
        elseif nargin==4
            obj.X = X;
            obj.genelist=genelist;
            obj.s=s;
            obj.c=c;
        end        
    end
    
      function r = libsz(obj)
         r = sum([obj.X]);
      end
      function r = norm(obj)
         r = sc_norm(obj.X);
      end
      
      function c_cell_idx=get.c_cell_idx(obj)
          c_cell_idx=1:size(obj.X,2);
      end
   function obj = set.c_cell_idx(obj,~)
      fprintf('%s%d\n','c_cell_idx is: ',obj.c_cell_idx(1:5))
      error('You cannot set the c_cell_idx property');
   end      

   function disp(td)
      fprintf(1,...
         '%d genes x %d cells\n',...
         size(td.X,1),size(td.X,2));
   end
   
   end
% https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html   
end


classdef SingleCellExperiment
   properties
      X double {mustBeNumeric, mustBeFinite}
      genelist 
      drmethod {mustBeMember(drmethod,{'tsne','umap','phate'})} = 'tsne'
      s (:,3)
      c
   end
   methods
    function obj = SingleCellExperiment(val,b)
        if nargin == 1
            obj.X = val;
        elseif nargin==2
            obj.X = val;
            obj.genelist=b;
        end
        
    end
    
      function r = libsz(obj)
         r = sum([obj.X]);
      end
      function r = norm(obj)
         r = sc_norm(obj.X);
      end
   end
end
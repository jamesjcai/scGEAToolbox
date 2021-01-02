classdef SingleCellEmbedding
   properties     
      s                               % 
      name {mustBeMember(name,{'tSNE','UMAP','PHATE','Unnamed'})} = 'Unnamed'
      method EmbeddingMethods
   end   
   properties (Dependent)
      NumDimensions
   end 

   methods
   function obj = SingleCellEmbedding(s,name)
        if nargin<2, name="Unnamed"; end
        obj.s = s;
        obj.name=name;
   end
   
   function m = get.NumDimensions(obj)   
      m = size(obj.s,2); 
   end
   
   function obj = set.NumDimensions(obj,~)
      fprintf('%s%d\n','NumDimensions is: ',obj.NumDimensions)
      error('You cannot set NumGenes property'); 
   end
 
   function p=plot(obj)
        p=gui.i_gscatter3(obj.s);
        title(obj.name);
   end
   end
end

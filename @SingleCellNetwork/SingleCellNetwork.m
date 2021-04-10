classdef SingleCellNetwork
   properties
      G                % digraph
   end
   
   properties (Dependent)
      NumGenes
      A
      g
   end 

   methods
   
   s=katz_score(obj,A,b)
   s=katz_centrality(obj,a,b)
   
   function obj = SingleCellNetwork(A,g)
        if nargin<1, A=[]; end
        if nargin<2 || isempty(g), g=string(transpose(1:size(A,1))); end
        g = matlab.lang.makeUniqueStrings(g);
        obj.G = digraph(A,g,'omitselfloops');
   end
   
   function A = get.A(obj)
       A=adjacency(obj.G,'weighted');
   end
   
   function g = get.g(obj)
       g = string(obj.G.Nodes.Variables);
   end
   
   function m = get.NumGenes(obj)   
      m = size(obj.G.Nodes,1); 
   end

   function obj = set.NumGenes(obj,~)
      fprintf('%s%d\n','NumGenes is: ',obj.NumGenes)
      error('You cannot set NumGenes property'); 
   end
 
    function p = plot(obj)
        p=plot(obj.G);
    end
    
    function p = plotweighted(obj)
        p=plot(obj.G);
        obj.G.Edges.LWidths = 7*obj.G.Edges.Weight/max(abs(obj.G.Edges.Weight));
        p.LineWidth = abs(obj.G.Edges.LWidths);        
    end
    
    function obj = makesparse(obj,q)
        if nargin<2, q=0.95; end
        a=max(abs(obj.A(:)));
        if a==0, a=1; end
        obj.A=obj.A./a;
        obj.A=obj.A.*(abs(obj.A)>quantile(abs(obj.A(:)),q));
        if ~issparse(obj.A)
            obj.A=sparse(obj.A);
        end
    end  
   end  

end

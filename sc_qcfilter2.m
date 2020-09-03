function [X,genelist,keptidxv]=sc_qcfilter2(X,genelist,varargin)

    if nargin<2
        error(sprintf('USAGE: [x,g]=sc_qcfilter2(X,genelist);'));
    end
    % if nargin<2, genelist=string(num2cell(1:size(X,1)))'; end
    
   p = inputParser;
   addOptional(p,'removemtgenes',true,@islogical);
   addOptional(p,'smplmethod',"Jackknife",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["jackknife","bootstrap"]));
   addOptional(p,'tdmethod',"CP",@(x) (isstring(x)|ischar(x))&ismember(upper(string(x)),["CP","TUCKER"]));
   addOptional(p,'min_cells_nonzero',0.01,@(x) isnumeric(x) & x>0);
   addOptional(p,'dropout',0.01,@(x) isnumeric(x) & x>0 & x<1);
   addOptional(p,'mtratio',0.1,@(x) isnumeric(x) & x>0 & x<1);
   addOptional(p,'libsize',1000,@(x) isnumeric(x) & x>10);
  
   parse(p,varargin{:});
   libsize=p.Results.libsize;
   mtratio=p.Results.mtratio;
   dropout=p.Results.dropout; 
   min_cells_nonzero=p.Results.min_cells_nonzero;
   removemtgenes=p.Results.removemtgenes;
   
[X,keptidx]=sc_rmmtcells(X,genelist,mtratio);
keptidxv{1}=keptidx;
if removemtgenes
    [X,genelist]=sc_rmmtgenes(X,genelist);
end
oldsz=0;
newsz=1;
c=1;
while ~isequal(oldsz,newsz)
    oldsz=size(X);
    [X,genelist]=sc_filterg(X,genelist,dropout);
    [X,keptidx]=sc_filterc(X);
    keptidxv{end+1}=keptidx;
    [X,genelist]=sc_selectg(X,genelist,1,min_cells_nonzero);
    [X,keptidx]=sc_selectc(X,libsize);
    keptidxv{end+1}=keptidx;
    newsz=size(X);
    c=c+1;
end

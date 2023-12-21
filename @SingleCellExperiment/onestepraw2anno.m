function obj = onestepraw2anno(obj, speciesid)

if nargin < 2 || isempty(speciesid), speciesid = 'human'; end
obj = obj.qcfilter;        
obj = obj.embedcells('tsne3d',true);        
obj = obj.clustercells([], [], true);
obj = obj.assigncelltype(speciesid, false);
end


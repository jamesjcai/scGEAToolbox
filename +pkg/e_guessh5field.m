function [y]=e_guessh5field(filnm,prefx,varm,required)

if nargin<4, required=false; end

y=[];
if iscell(varm) || ~isstring(varm), varm=string(varm); end
if iscell(prefx) || ~isstring(prefx), prefx=string(prefx); end

for l=1:length(prefx)
for k=1:length(varm)
    try
        y=h5read(filnm,sprintf('%s%s',prefx(l),varm(k)));
    catch
        continue;
    end
    break;
end
end

if isempty(y) && required
    %errorStruct.message = 'Data file not found.';
    %errorStruct.identifier = 'MyFunction:fileNotFound';
   error('Required value is missing.');
end
end

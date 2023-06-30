function [y]=e_guessh5field(filnm,prefx,varm)

y=[];
if iscell(varm), varm=string(varm); end
for k=1:length(varm)
    try
        y=h5read(filnm,sprintf('%s%s',prefx,varm(k)));
    catch
        continue;
    end
    break;
end
end

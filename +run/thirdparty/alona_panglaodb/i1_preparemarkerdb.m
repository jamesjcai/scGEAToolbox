T=readtable('PanglaoDB_markers_27_Mar_2020.tsv','filetype','text');

sp="Mm"; 
%sp="Hs";

i=contains(string(T.species),sp);
T=T(i,:);
organlist=string(T.organ);

%%
organ="All";
% organ="Heart";
% organ="Pancreas";
% organ="Brain";
% organ="Immune system";

% organ="Thymus"

i=organ==organlist;
if isempty(i), error('No organ.'); end


if strcmpi(organ,"All")
    disp('All organs.');
    outfile=sprintf('markerlist_%s_panglaodb.txt',lower(sp));
else   % if strcmpi(organ,"Heart") && ~isempty(i)    
    T=T(i,:);
    organv=lower(matlab.lang.makeValidName(organ));
    if ~exist(organv,'dir'), mkdir(organv); end
    outfile=sprintf('%s/markerlist_%s_panglaodb.txt',organv,lower(sp));
end



a=string(unique(T.cellType));
gt=string(T.cellType);
genesymbollist=string(T.officialGeneSymbol);
organlist=string(T.organ);

fid=fopen(outfile,'w');
for k=1:length(a)
    k;
    idx=find(gt==a(k));
    fprintf(fid,'%s\t',a(k));
    % fprintf(fid,'%s\t',organlist(k));
    fprintf(fid,'%s,', genesymbollist(idx));
    fprintf(fid,'\n');
end
fclose(fid);

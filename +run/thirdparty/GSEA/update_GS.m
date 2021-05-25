function update_GS(data_name)
% Function for reading GeneSet databaset from .xls file and storing as .mat
% file.
% Input:
% data_name - filename of GS database with filename extension (.xls,.xlsx, etc.)

dot_ind = find(data_name=='.');
if isempty(dot_ind)
    error('No file extension included.')
end

%load data
[tmp1,~,tmp2] = xlsread(data_name);
disp([data_name(1:dot_ind-1) ' dataset loaded.'])

GS.ID = tmp2(:,1);        %GeneSet ID
GS.descr = tmp2(:,2);       %GeneSet descriptions
GS.nb = length(GS.ID);    %Number of GS in databaset

GS.entrez = cell(GS.nb,1);      %Gene entrez ID within each GS
GS.entrez_nb = zeros(GS.nb,1);  %Number of genes within each GS
for a=1:GS.nb
    disp(['GS no. ' num2str(a) '/' num2str(GS.nb)])
    tmp = tmp1(a,3:end);
    GS.entrez{a} = unique(tmp(~isnan(tmp)));
    GS.entrez_nb(a) = length(GS.entrez{a});
    
    if isnumeric(GS.ID{a})
        GS.ID{a} = num2str(GS.ID{a});
    end        
end

save([data_name(1:dot_ind-1) '.mat'],'GS')
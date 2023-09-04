function [output]=i_symbol2ncbiid(input)
if nargin<1
    input=["AKAIN1","AKAP1"];
end
output=zeros(length(input),1);
pw1=fileparts(mfilename('fullpath'));
idfile=fullfile(pw1,'..','resources','genename2ncbigeneid.xlsx');
% if ~exist(idfile,'file')
%     options = weboptions('Timeout',21);    
%     disp('Downloading genename2ncbigeneid.xlsx...');
%     websave(idfile,'url',options);    
% end
T=readtable(idfile);
[y,idx]=ismember(input,string(T.ApprovedSymbol));
output(y)=T.NCBIGeneID(idx(y));
end

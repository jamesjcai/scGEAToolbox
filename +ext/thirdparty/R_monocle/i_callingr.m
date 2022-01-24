
% csvwrite('in.csv',[originaldat]);
if ~exist('X.csv','file')
    csvwrite('X.csv',X);
end
if ~exist('genelist.txt','file')
    writetable(table(genelist),'genelist.txt');
end

Rpath = 'C:\Program Files\R\R-3.6.0\bin';
RscriptFileName = 'Z:\Cailab\mouse_neurons\adult_P10_cortex_SRR6061129\monocleMatlab.R';
RunRcode(RscriptFileName)

% !R CMD BATCH monocleMatlab.R
dat = csvread('pt_X.csv',1,1);


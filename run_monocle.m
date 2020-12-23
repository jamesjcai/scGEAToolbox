function [t,s]=run_monocle(X,plotit)
% Run Monocle pseudotime analysis

%[t_mono,s_mono]=run_monocle(X,true);

if isempty(FindRpath)
   error('Rscript.ext is not found.');
end
if nargin<2, plotit=false; end
oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/R_monocle');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('output.csv','file')
    delete('output.csv');
end
if issparse(X), X=full(X); end
csvwrite('input.csv',X);
% Rpath = 'C:\Program Files\R\R-3.6.0\bin';
% RscriptFileName = 'Z:\Cailab\mouse_neurons\adult_P10_cortex_SRR6061129\monocleMatlab.R';
% RunRcode('monocleMatlab_3d.R');
RunRcode('script.R');
if exist('output.csv','file')
    dat = csvread('output.csv',1,1);
    t=dat(:,1);
    s=dat(:,2:end);
else
    t=[];
    s=[];
end
if exist('input.csv','file')
    delete('input.csv');
end
if exist('output.csv','file')
    delete('output.csv');
end
cd(oldpth);

if plotit
    if ~isempty(s) && ~isempty(t)
        % i_myscatter(s,t);
        i_gscatter3(s,t);
    end
end




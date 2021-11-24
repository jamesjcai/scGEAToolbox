function [t,s]=monocle(X)
% Run Monocle pseudotime analysis
%[t_mono,s_mono]=run.monocle(X,true);

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_monocle');
if ~isok, error(msg); t=[]; s=[]; return; end

if ~isdebug
if exist('./imput.mat','file'), delete('./imput.mat'); end
if exist('./output.csv','file'), delete('./output.csv'); end
end

save('input.mat','X','-v7.3');
% writematrix(X,'input.csv');
% Rpath = 'C:\Program Files\R\R-3.6.0\bin';
% RscriptFileName = 'Z:\Cailab\mouse_neurons\adult_P10_cortex_SRR6061129\monocleMatlab.R';
% pkg.RunRcode('monocleMatlab_3d.R');
pkg.RunRcode('script.R');
if exist('./output.csv','file')
    dat=readmatrix('output.csv');
    t=dat(:,2);
    s=dat(:,3:end);
else
    t=[];
    s=[];
end


if ~isdebug
if exist('./imput.mat','file'), delete('./imput.mat'); end
if exist('./output.csv','file'), delete('./output.csv'); end
end
cd(oldpth);

% if plotit
%     if ~isempty(s) && ~isempty(t)
%         % i_myscatter(s,t);
%         i_gscatter3(s,t);
%     end
% end

end
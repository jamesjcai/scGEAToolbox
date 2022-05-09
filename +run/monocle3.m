function [t,s]=monocle3(X)
% Run Monocle3 pseudotime analysis
%[t_mono,s_mono]=run.monocle3(X);

    t=[]; s=[];
    isdebug=false;
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_monocle3');
    if ~isok, error(msg); end
    
    tmpfilelist={'input.h5','output.h5'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    pkg.e_writeh5(X,[],'input.h5');
    pkg.RunRcode('script.R');
    if exist('./output.h5','file')
        t = h5read('output.h5','/t');
        s = h5read('output.h5','/s');
    %    dat=readmatrix('output.csv');
    %    t=dat(:,2);
    %    s=dat(:,3:end);
    end
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
end

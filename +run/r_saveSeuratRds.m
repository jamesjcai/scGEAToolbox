function [status]=r_saveSeuratRds(sce,filename)
    [status]=0;
    isdebug=false;
    if nargin<2, error('run.saveSeuratRds(sce,filename)'); end
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_SeuratSaveRds');
    if ~isok, error(msg); return; end

    tmpfilelist={'input.h5','output.Rds'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    pkg.e_writeh5(sce.X,sce.g,'input.h5');
    %sc_writefile('input.txt',sce.X,sce.g);
    %    if isdebug, return; end
    Rpath=getpref('scgeatoolbox','rexecutablepath');
    pkg.RunRcode('script.R',Rpath);
    
    [status]=copyfile('output.Rds',filename,'f');
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
end
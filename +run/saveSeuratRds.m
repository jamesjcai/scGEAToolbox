function [status]=saveSeuratRds(sce,filename)
    [status]=0;
    if nargin<2, error('run.saveSeuratRds(sce,filename)'); end
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_SeuratSaveRds');
    if ~isok, error(msg); end
    if exist('output.Rds','file'), delete('output.Rds'); end
    sc_writefile('input.txt',sce.X,sce.g);
    RunRcode('script.R');
    [status]=copyfile('output.Rds',filename,'f');
    if exist('input.txt','file'), delete('input.txt'); end
    if exist('output.Rds','file'), delete('output.Rds'); end
    cd(oldpth);
end

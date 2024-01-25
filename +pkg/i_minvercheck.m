function [needupdate]=i_minvercheck
% min version update check
needupdate = false;
olddir = pwd();
cdgea;
try
    Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox');
    if ~exist('TIMESTAMP','file')
        fid=fopen('TIMESTAMP','w');
        fprintf(fid,'%s\n',Col.pushed_at);
        fclose(fid);
    else
        fid=fopen('TIMESTAMP','r');
        d=textscan(fid,'%s');
        d=d{1}{1};
        fclose(fid);
        if ~strcmp(d,Col.pushed_at)

            answer=questdlg('A minor update is now available. Learn how to install it?','', ...
                'Yes','Remind me later','Skip this update','Skip this update');
            switch answer
                case 'Yes'
                    web('https://scgeatoolbox.readthedocs.io/en/latest/quick_installation.html');
                    needupdate=true;
                    return;
                case 'Remind me later'
                case 'Skip this update'
                    fid=fopen('TIMESTAMP','w');
                    fprintf(fid,'%s\n',Col.pushed_at);
                    fclose(fid);
                otherwise
            end        
        end
    end
catch
end
cd(olddir);
end

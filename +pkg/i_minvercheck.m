function [needupdate] = i_minvercheck(parentfig)
if nargin<1, parentfig = []; end
% min version update check
needupdate = false;
mfolder = fileparts(mfilename('fullpath'));

try
    Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox');

    stfile = fullfile(mfolder,'..','TIMESTAMP');

    if ~exist(stfile,'file')
        fid=fopen(stfile,'w');
        fprintf(fid,'%s\n',Col.pushed_at);
        fclose(fid);
        needupdate = true;        
    else
        fid=fopen(stfile,'r');
        d=textscan(fid,'%s');
        d=d{1}{1};
        fclose(fid);
        if ~strcmp(d,Col.pushed_at)

            answer=gui.myQuestdlg(parentfig, 'A minor update is now available. Learn how to install it?','', ...
                {'Yes','Remind me later','Skip this update'},'Skip this update');
            switch answer
                case 'Yes'
                    gui.gui_uishowrefinfo('Quick Installation', parentfig);
                    needupdate = true;
                    return;
                case 'Remind me later'
                case 'Skip this update'
                    fid=fopen(stfile,'w');
                    fprintf(fid,'%s\n',Col.pushed_at);
                    fclose(fid);
                otherwise
            end
        else
            gui.myHelpdlg(parentfig, 'No update is available.','');
        end
    end
catch ME
    errordlg(ME.message);
end

end

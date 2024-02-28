function [needupdate]=i_minvercheck
% min version update check
needupdate = false;
mfolder = fileparts(mfilename('fullpath'));

try
    Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox');

    stfile =   fullfile(mfolder,'..','TIMESTAMP');

    if ~exist(stfile,'file')
        fid=fopen(stfile,'w');
        fprintf(fid,'%s\n',Col.pushed_at);
        fclose(fid);
    else
        fid=fopen(stfile,'r');
        d=textscan(fid,'%s');
        d=d{1}{1};
        fclose(fid);
        if ~strcmp(d,Col.pushed_at)

            answer=questdlg('A minor update is now available. Learn how to install it?','', ...
                'Yes','Remind me later','Skip this update','Skip this update');
            switch answer
                case 'Yes'
                    %web('https://scgeatoolbox.readthedocs.io/en/latest/quick_installation.html');

                    prompt = {'Copy the following code and run it in MATLAB:'};
                    dlgtitle = 'Quick Installation/Upgrade';
                    fieldsize = [18 75];
                    definput = {sprintf('tic;\ndisp(''Installing scGEAToolbox...'')\nunzip(''https://github.com/jamesjcai/scGEAToolbox/archive/main.zip'');\naddpath(''./scGEAToolbox-main'');\ntoc;\nif exist(''scgeatool.m'',''file'')\n    disp(''scGEAToolbox installed!'')\nend\nsavepath(fullfile(userpath,''pathdef.m''));')};
                    inputdlg(prompt,dlgtitle,fieldsize,definput);

                    needupdate=true;
                    return;
                case 'Remind me later'
                case 'Skip this update'
                    fid=fopen(stfile,'w');
                    fprintf(fid,'%s\n',Col.pushed_at);
                    fclose(fid);
                otherwise
            end
        else
            waitfor(helpdlg('No update is available.',''));
        end
    end
catch ME
    errordlg(ME.message);
end

end

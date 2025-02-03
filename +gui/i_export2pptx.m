function i_export2pptx(F, glist)

if nargin < 2, glist = {[]}; end

if ~license('test','matlab_report_gen') && ~isempty(which('mlreportgen.report.Report'))
    warndlg('Unable to check out a Report Generator license. This function requires MATLAB Report Generator.', '');
    return;
end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'resources', 'Misc', 'myTemplate.pptx');
%dbfile1 = fullfile(pw1, '+run', 'external', 'stringdb', 'stringdb_human.mat');

import mlreportgen.ppt.*;
answer = questdlg('Export to PowerPoint?');
switch answer
    case 'Yes'
        if ~usejava('desktop')
           warndlg('This function is not supported in MATLAB Online and may not work properly in standalone applications.', '');
           return;
        end        
        if ismcc || isdeployed
            makePPTCompilable();
            warndlg('This function may not work properly in standalone applications.', '');
        end
        fw = gui.gui_waitbar;
        N = length(F);
        images = cell(N, 1);

        OUTppt = [tempname, '.pptx'];
        ppt = Presentation(OUTppt, pth);

%         try
%     import mlreportgen.report.*; % Example usage of Report Generator
%     rpt = Report('MyReport', 'pdf');
%     fprintf('MATLAB Report Generator is working correctly.\n');
% catch ME
%     if contains(ME.message, 'MATLAB Report Generator is not installed')
%         warning('MATLAB Report Generator is not installed or licensed.');
%     else
%         rethrow(ME); % If it's another error, rethrow it
%     end
% end

        try
            open(ppt);
        catch ME
            pause(0.5);
            gui.gui_waitbar(fw, true);
            pause(0.5);
            waitfor(errordlg(ME.message,''));
            return;
        end
        
        warning off
        for k = 1:N
            if isvalid(F{k})
                images{k} = [tempname, '.png'];
                
                saveas(F{k}, images{k});
                %images{k} = [tempname,'.emf'];
                %saveas(F{k},images{k},'meta');

                if ~isempty(glist{1})
                    slide3 = add(ppt, 'Small Title and Content');
                    replace(slide3, 'Title', glist(k));
                else
                    slide3 = add(ppt, 'Content Only');
                end
                replace(slide3, 'Content', Picture(images{k}));
            end
        end
        warning on
        close(ppt);
        len = length(images);
        for i = 1:len
            delete(images{i});
        end
        gui.gui_waitbar(fw);
        %try
           rptview(ppt);
        %catch ME
        %    warndlg(ME.message,'');
        %end
    otherwise
       
end

end

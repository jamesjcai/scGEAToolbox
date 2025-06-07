function i_export2pptx(F, glist, parentfig)

if nargin<3, parentfig = []; end

if nargin < 2, glist = {[]}; end

if ~license('test','matlab_report_gen') && ~isempty(which('mlreportgen.report.Report'))
    gui.myWarndlg(parentfig, ['Unable to check out a Report ' ...
        'Generator license. This function requires MATLAB ' ...
        'Report Generator.']);
    return;
end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');
%dbfile1 = fullfile(pw1, '+run', 'external', 'stringdb', 'stringdb_human.mat');

import mlreportgen.ppt.*;
answer = gui.myQuestdlg(parentfig, 'Export to PowerPoint?');
switch answer
    case 'Yes'
        if ~usejava('desktop')
           gui.myWarndlg(parentfig, ['This function is not supported in ' ...
               'MATLAB Online and may not work properly in ' ...
               'standalone applications.']);
           return;
        end        
        if ismcc || isdeployed
            makePPTCompilable();
            gui.myWarndlg(parentfig, ['This function may not work properly in ' ...
                            'standalone applications.']);
        end
        fw = gui.myWaitbar(parentfig);
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
            gui.myWaitbar(parentfig, fw, true);
            pause(0.5);
            gui.myErrordlg(parentfig, ME.message);
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
        gui.myWaitbar(parentfig, fw);
        rptview(ppt);
    otherwise
       
end

end


%{

import mlreportgen.ppt.*

pptFile = 'myPresentation.pptx';  % PowerPoint file name

% Check if PowerPoint file exists, else create a new one
if isfile(pptFile)
    ppt = Presentation(pptFile); % Open existing PPTX
else
    ppt = Presentation(pptFile); % Create new PPTX
end

% Loop to create and add multiple figures
for k = 1:3
    % Create a figure (replace with your plotting code)
    fig = figure;
    plot(rand(10,1), 'LineWidth', 2);
    title(['Vector Figure ' num2str(k)]);
    
    % Save figure as EMF (vector format)
    emfFile = ['figure' num2str(k) '.emf'];
    print(fig, emfFile, '-dmeta'); % Save as EMF
    
    % Add slide to PowerPoint
    slide = add(ppt, 'Title and Content');
    replace(slide, 'Title', ['Vector Figure ' num2str(k)]);
    add(slide, Picture(emfFile));

    % Close the figure
    close(fig);
end

% Save and close the presentation
close(ppt);

disp('Figures added to PowerPoint successfully with EMF format.');
%}
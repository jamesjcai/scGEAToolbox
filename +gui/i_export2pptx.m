function i_export2pptx(F, glist, parentfig)
%I_EXPORT2PPTX Put figures into a PowerPoint deck, one slide each.
%
%   gui.i_export2pptx(F, glist, parentfig)
%
%   F      cell array of figure handles, one slide per figure
%   glist  matching slide titles, {''} for untitled slides
%
%   A figure holding a tab group is offered a slide per tab instead of a
%   single slide. GETFRAME photographs what is on screen, so tabs that were
%   not selected used to be absent from the deck entirely although they were
%   fully drawn: a gui.MYFIGURE from one of the SC_UITABGRPFIG_* plotters
%   carries a tab per gene, and only the gene on top was ever exported.
%   GUI.I_CAPTUREFRAME draws each tab on its own, so a slide is the plots
%   rather than a photograph of the window around them.
%
%   See also gui.i_savemainfig, gui.i_save2pptx, gui.sc_uitabgrpfig_expplot.

if nargin<3, parentfig = []; end
if nargin < 2, glist = {[]}; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

[hasReportGen, msg] = pkg.i_isreportgenavailable('ppt');
if ~hasReportGen
    gui.myWarndlg(parentfig, sprintf('%s This function requires MATLAB Report Generator.', msg));
    return;
end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');

import mlreportgen.ppt.*;

% What to photograph, and under what title. One entry per slide; the question
% about tabs is asked in there, so nothing below has to care whether this is
% one figure or six tabs of one.
[capfig, captab, captitle] = in_planslides(F, glist, parentfig);
if ~isempty(capfig)
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
        % Show the deck instead of asking where to put it first. The
        % point of the button is to see the slides, and PowerPoint's own
        % Save As decides where they end up; a UIPUTFILE here also made
        % this button behave unlike GUI.I_SAVE2PPTX, which opens the deck
        % straight away.
        %
        % What made the original version unrecoverable was TEMPNAME: a
        % random basename directly in TEMP, so a deck that never launched
        % was a deck the user could not find, and neither RPTVIEW nor
        % WINOPEN raises when the shell declines to open a file.
        % PKG.I_TEMPDIRFILE gives a named per-process folder and a
        % timestamped file instead, which the fallback at the end can
        % point at.
        [~, OUTppt] = pkg.i_tempdirfile('scgeatool_pptx', 'pptx');

        fw = gui.myWaitbar(parentfig);
        N = numel(capfig);
        images = cell(N, 1);

        % Drawing a tab does not select it, so the loop below normally
        % leaves the figure exactly where the user had it. The one path that
        % does click - GUI.I_CAPTUREFRAME falling back to GETFRAME for a tab
        % EXPORTGRAPHICS will not draw - is why this is still here.
        restoretab = in_holdtabselection(captab);   %#ok<NASGU>

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

        % Capture and restore rather than a bare off/on pair: the pair leaves
        % warnings disabled for the rest of the session if anything between the
        % two lines throws, and its 'on' re-enables warnings the caller may have
        % silenced deliberately instead of restoring what they had.
        warnState = warning();
        restoreWarn = onCleanup(@() warning(warnState));
        warning('off', 'all');
        for k = 1:N
            if pkg.i_isvalid(capfig{k})
                images{k} = [tempname, '.png'];
                gui.i_captureframe(capfig{k}, captab{k}, images{k});

                % option 1
                % saveas(F{k}, images{k});

                % option 2
                % exportapp(F{k}, images{k});

                % option 3
                %
                % set(F{k}, 'Color', 'white');
                % tabGroups = findall(F{k}, 'Type', 'uitabgroup');
                % for i = 1:length(tabGroups)
                %     try
                %         set(tabGroups(i), 'BackgroundColor', 'white');
                %     catch
                %         disp("error: set(tabGroups(i), 'BackgroundColor', 'white');");
                %         try
                %             set(tabGroups(i), 'ForegroundColor', 'white');
                %         catch
                %             disp("error: set(tabGroups(i), 'ForegroundColor', 'white');");
                %         end
                %     end
                % end
                % tabs = findall(F{k}, 'Type', 'uitab');
                % for i = 1:length(tabs)
                %     try
                %         set(tabs(i), 'BackgroundColor', 'white');
                %     catch
                %     end
                % end

                % images{k} = [tempname,'.emf'];
                % saveas(F{k},images{k},'meta');

                if strlength(captitle(k)) > 0
                    slide3 = add(ppt, 'Small Title and Content');
                    replace(slide3, 'Title', char(captitle(k)));
                else
                    slide3 = add(ppt, 'Content Only');
                end
                replace(slide3, 'Content', Picture(images{k}));
            end
        end
        close(ppt);
        len = length(images);
        for i = 1:len
            delete(images{i});
        end

        % Put the caller's warning state back before opening the file.
        % RESTOREWARN above does it at function exit, which is after this
        % point, so anything PKG.I_OPENOUTPUTFILE had to say about failing
        % to launch the deck went into a session where warnings were still
        % off, and was never seen.
        clear restoreWarn;

        gui.myWaitbar(parentfig, fw);
        if ~pkg.i_openoutputfile(OUTppt)
            gui.myHelpdlg(parentfig, sprintf( ...
                ['The presentation could not be opened here. ', ...
                'It has been saved as\n\n%s'], OUTppt));
        end
end

end


function [capfig, captab, captitle] = in_planslides(F, glist, parentfig)
% Ask what the deck should contain, then hand the answer to
% GUI.I_PPTXSLIDEPLAN. An empty CAPFIG means the user declined and nothing
% should be written.

capfig = {};
captab = {};
captitle = strings(0, 1);

tabs = gobjects(0);
if isscalar(F), tabs = gui.i_figtabs(F{1}); end

if numel(tabs) > 1
    % Asked, not assumed: a plotter handed twenty genes makes twenty tabs,
    % and twenty slides is not always the wanted answer.
    answer = gui.myQuestdlg(parentfig, sprintf(['This figure has %d tabs. ' ...
        'Export every tab as its own slide, or only the tab on top?'], ...
        numel(tabs)), 'Export to PowerPoint', ...
        {'All tabs', 'Current tab only', 'Cancel'}, 'All tabs');
    switch answer
        case 'All tabs'
            [capfig, captab, captitle] = gui.i_pptxslideplan(F, glist, tabs);
        case 'Current tab only'
            [capfig, captab, captitle] = gui.i_pptxslideplan(F(1), {''});
        otherwise
            % Cancelled, or the dialog was dismissed.
    end
    return;
end

if ~strcmp(gui.myQuestdlg(parentfig, 'Export to PowerPoint?'), 'Yes')
    return;
end
[capfig, captab, captitle] = gui.i_pptxslideplan(F, glist);
end


function [c] = in_holdtabselection(captab)
% An ONCLEANUP that puts the selected tab back, or an inert one when this
% export is not touching tabs at all.

c = onCleanup(@() []);
hit = find(~cellfun(@isempty, captab), 1);
if isempty(hit), return; end

tg = captab{hit}.Parent;
if ~pkg.i_isvalid(tg), return; end

wastab = tg.SelectedTab;
c = onCleanup(@() in_reselect(tg, wastab));
end


function in_reselect(tg, tab)
if pkg.i_isvalid(tg) && pkg.i_isvalid(tab)
    tg.SelectedTab = tab;
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

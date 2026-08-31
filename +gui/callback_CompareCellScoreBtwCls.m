function [needupdatesce] = callback_CompareCellScoreBtwCls(src, ~)

%{
1. Input:
   - `sce`: An object containing the cellular data.
   - `glist`: A list of gene/region names (optional).
   - `thisc`: The current cell group being analyzed.
   - `selecteditem`: Determines which type of score to calculate or display.

2. Output:
    - Displays a graph showing the calculated scores along with labels for each dataset.
%}

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

needupdatesce = false;

% Tracks which of the original cells are still in play. Selecting two groups
% below subsets sce, so scores would then cover only those cells; keepidx
% lets the save step recognize that case and skip saving.
keepidx = true(sce.NumCells, 1);

aa = 'Yes, compare scores (violinplot)';
bb = 'No, just show values';

answer2 = gui.myQuestdlg(FigureHandle, sprintf(['This function will calculates a score for each cell. After the scores are calculated, do you want to ', ...
    'compare score values between different cell groups?']), '', ...
    {bb, aa}, bb);
switch answer2
    case aa
        showcomparision = true;
    case bb
        showcomparision = false;
    otherwise
        return;
end

if ~showcomparision
    thisc = ones(sce.NumCells, 1);
else
    allowunique = false;
    % Several grouping variables may be picked here; they are
    % crossed into one composite label per cell ("Macrophages |
    % IL"). Everything downstream treats thisc as a per-cell
    % label vector, so the composite needs no special handling
    % past this point.
    [thisc, clabel] = gui.i_selectnclass(sce, allowunique, ...
        [], [], FigureHandle);
    if isempty(thisc), return; end
    if isscalar(unique(thisc))
        answer = gui.myQuestdlg(FigureHandle, "All cells are in the same group. No comparison will be made. Continue?", ...
            "", {'Yes', 'No', 'Cancel'}, 'No', 'warning');
        switch answer
            case 'Yes'
            otherwise
                return;
        end
    else    % length(unique(thisc)) ~= 1
        % Shared with gui.callback_Violinplot so both group choosers
        % look and behave the same - see gui.i_selectgroupsubset.
        idx2 = gui.i_selectgroupsubset(thisc, clabel, FigureHandle);
        if isempty(idx2), return; end
        sce = sce.selectcells(idx2);  % OK
        thisc = thisc(idx2);
        keepidx = idx2(:);
    end
end
drawnow;
[selecteditem, speciestag] = gui.i_selgenesetcollection(FigureHandle);
if isempty(selecteditem), return; end

posg = []; ttxt = [];
% Which scoring algorithm ends up being used. Some
% branches pick it here, others let gui.e_cellscore ask;
% either way it has to reach the figure, or the plot
% cannot say how its numbers were produced.
methodid = [];
switch selecteditem
    % case 'Global Coordination Level (GCL) [PMID:33139959]'

    case 'Define a New Score...'
        ttxt = 'Customized Score';

        if gui.i_isuifig(FigureHandle)
            newcstype = gui.myInputdlg({'Name of the new cell score:'}, ...
                '', {ttxt}, FigureHandle);
        else
            newcstype = inputdlg('Name of the new cell score:', ...
                '', [1, 50], {ttxt});
        end

        if isempty(newcstype), return; end
        ttxt = newcstype;
        [posg] = gui.i_selectngenes(sce.g, [], FigureHandle);
        if isempty(posg)
            gui.myHelpdlg(FigureHandle, 'No feature genes selected.', '')
            return;
        end
        [y, methodid] = gui.e_cellscore(sce, posg, [], [], FigureHandle);
        % Saving is handled centrally after the switch, so every
        % score source offers it on the same terms.
    case 'MSigDB Molecular Signatures'
        % speciestag = gui.i_selectspecies(2, true, FigureHandle);
        % if isempty(speciestag), return; end
        try
            [posg, ctselected] = gui.i_selectMSigDBGeneSets(speciestag, false, FigureHandle);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        ttxt = ctselected;
        if isempty(posg) || isempty(ctselected), return; end

        n = length(posg);
        y = cell(n,1);
        [~, methodid] = gui.i_pickscoremethod([], FigureHandle);
        if isempty(methodid), return; end

        fw = gui.myWaitbar(FigureHandle);
        for k = 1:n
            y{k} = gui.e_cellscore(sce, posg{k}, ...
                methodid, false, FigureHandle);
        end
        gui.myWaitbar(FigureHandle, fw);


    case 'PanglaoDB Cell Type Markers'
        if isempty(speciestag)
            speciestag = gui.i_selectspecies(2, true, FigureHandle);
        end
        if isempty(speciestag), return; end

        pw1 = fileparts(mfilename('fullpath'));
        pth = fullfile(pw1, '..', 'external', 'fun_alona_panglaodb');
        switch lower(speciestag)
            case {'human', 'hs'}
                speciestag_short = 'hs';
            case {'mouse', 'mm'}
                speciestag_short = 'mm';
            otherwise
                speciestag_short = 'hs';
        end

        markerfile = fullfile(pth, sprintf('marker_%s.mat', speciestag_short));
        if exist(markerfile, 'file')
            load(markerfile, 'Tm');
        else
            Tm = readtable(fullfile(pth, sprintf('markerlist_%s.txt', speciestag_short)), ...
                'ReadVariableNames', false, 'Delimiter', '\t');
        end

        ctlist = string(Tm.Var1);
        listitems = sort(ctlist);
        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, listitems, ...
                'Select Class:');
        else
            [indx, tf] = listdlg('PromptString', ...
                {'Select Class'}, ...
                'SelectionMode', 'single', 'ListString', listitems, ...
                'ListSize', [220, 300]);
        end
        if ~tf == 1, return; end
        ctselected = listitems(indx);
        % idx=find(matches(ctlist,ctselected));
        idx = matches(ctlist, ctselected);
        ctmarkers = Tm.Var2{idx};
        posg = string(strsplit(ctmarkers, ','));
        posg(strlength(posg) == 0) = [];
        ttxt = ctselected;
        if isempty(posg) || isempty(ctselected), return; end

        [y, methodid] = gui.e_cellscore(sce, posg, [], [], FigureHandle);

    case {'Predefined Custom Gene Sets', 'Glycobiology Gene Sets'}
        % if ~gui.gui_showrefinfo('Predefined Cell Score', FigureHandle), return; end
        % Both entries score the same collection: pkg.e_cellscores already
        % appends pkg.e_glycogenesets to its table, so the glycobiology
        % entry is that table filtered to SignatureTag='Glycobiology'.
        [y, ttxt, posg, methodid] = i_scorePredefinedSets(sce, ...
            selecteditem, FigureHandle);
        % TF activity analysis has moved to Analyze > TF Activity Analysis
        % (gui.callback_TFActivity), which provides a dedicated workflow
        % with TF pre-selection and multi-TF display support.

    otherwise
        return;
end
% catch ME
%    %gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
%    return;
% end

if isempty(y), return; end

% Offer to keep the computed score(s) on the SCE. Declining only
% skips the save - the plots below are drawn either way. Nothing
% is offered when only a subset of cells was analysed.
[sce, needupdatesce] = i_saveScoresAsAttributes(sce, ...
    keepidx, y, ttxt, FigureHandle);

% assignin('base','y',y);
% assignin('base','thisc',thisc);

if showcomparision
    % if iscell(y)
    plotkind = i_pickComparisonPlotType(FigureHandle);
    if ~isempty(plotkind)
        % Say what the axes are and where the numbers came
        % from. The tab only names the score; without these
        % the reader cannot tell a UCell score from an
        % AddModuleScore one, or what the x groups are.
        labelinfo = struct();
        labelinfo.ylabel = 'Cell score';
        if exist('clabel', 'var') && ~isempty(clabel)
            labelinfo.xlabel = clabel;
        end
        labelinfo.method = i_methodlabel(methodid);
        gui.sc_uitabgrpfig_vioplot(y, ttxt, thisc, ...
            FigureHandle, plotkind, labelinfo);
    end
    % else
    %                    gui.i_violinplot(y, thisc, ttxt, true, [], posg, FigureHandle);
    %                    xlabel('Cell group');
    %                    ylabel('Cellular score');
    %                end
else
    %     [methodid]=gui.i_pickscatterstem('Scatter+Stem');
    %     if isempty(methodid), return; end
    %         f=gui.i_cascadefig(sce,glist(k),axx,bxx,k,methodid);
    %     [h1]=sc_scattermarker(sce.X,sce.g,...
    %                  sce.s,g,methodid);
    if iscell(y)
        % t=array2table(cell2mat({rand(10,1),rand(10,1),rand(10,1)}),'VariableNames',{'aa','bb','cc'});
        % assignin("base",'y',y);
        % assignin("base",'ttxt',ttxt);
        % assignin("base",'k',k);

        if length(y)>1
            % One tab per score for Heat/Stem; the score heatmap
            % is a different view (scores x cells, grouped).
            switch i_pickMultiScorePlotType(FigureHandle)
                case 'Heat'
                    gui.sc_uitabgrpfig_feaplot(y, string(ttxt), ...
                        sce.s, FigureHandle, 1);
                case 'Stem'
                    gui.sc_uitabgrpfig_feaplot(y, string(ttxt), ...
                        sce.s, FigureHandle, 2);
                case 'Score heatmap'
                    drawnow;
                    gui.i_scoreheatmap(cell2mat(y').', ttxt, ...
                        sce, FigureHandle, thisc);
            end
        elseif isscalar(y)
            figfun = i_pickScatterPlotType(FigureHandle);
            if ~isempty(figfun)
                figfun(sce, y{1}, posg, ttxt{1}, FigureHandle);
            end
        end
        % gui.sc_uitabgrpfig_expplot(y, markerlist, sce.s, FigureHandle, [axx, bxx]);
    else
        % Any saving already happened above, with the user asked.
        figfun = i_pickScatterPlotType(FigureHandle);
        if ~isempty(figfun)
            figfun(sce, y, posg, ttxt, FigureHandle);
        end
    end
end

if needupdatesce
    gui.myGuidata(FigureHandle, sce, src);
    if isa(src, 'matlab.apps.AppBase')
        src.sce = sce;
    end
end
end


function [y, ttxt, posg, methodid] = i_scorePredefinedSets(sce, ...
    collection, parentfig)
%I_SCOREPREDEFINEDSETS Pick signatures from the predefined table and score them.
%
%   collection decides which rows of the pkg.e_cellscores table are offered:
%     'Glycobiology Gene Sets'      - only the curated glycobiology modules
%     'Predefined Custom Gene Sets' - asks which collection(s) first, via
%                                     gui.i_picksignaturetags, then the
%                                     scores within them
%
%   These used to be two branches with two separate scoring calls, even
%   though pkg.e_cellscores appends pkg.e_glycogenesets to its own table and
%   tags those rows SignatureTag='Glycobiology' - the gene lists and the
%   resulting scores were identical either way. The glycobiology branch also
%   lacked this one's guard for signatures with too few expressed genes,
%   which left y and ttxt at different lengths and made
%   i_saveScoresAsAttributes silently skip saving every score, not just the
%   failed one. One code path with a row filter avoids both problems.
%
%   Returns an empty y when the user cancels or nothing could be scored; the
%   caller treats that as "stop here".

y = {}; ttxt = {}; posg = []; methodid = [];

[~, T] = pkg.e_cellscores(sce.X, sce.g, 0);

if strcmp(collection, 'Glycobiology Gene Sets')
    % Already a named collection - no point asking which one.
    rows = find(strcmp(string(T.SignatureTag), "Glycobiology"));
    prompt = 'Select Module:';
    dlgwidth = 300;
else
    [rows, taglabel] = gui.i_picksignaturetags(T, parentfig);
    if isempty(rows), return; end
    if strlength(taglabel) > 0
        prompt = char("Select Score (" + taglabel + "):");
    else
        prompt = 'Select Score:';
    end
    dlgwidth = 260;
end
if isempty(rows)
    gui.myWarndlg(parentfig, ...
        'No signatures are available in this collection.');
    return;
end
listitems = string(T.ScoreType(rows));

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, listitems, prompt);
else
    [indx, tf] = listdlg('PromptString', prompt, ...
        'SelectionMode', 'multiple', 'ListString', ...
        cellstr(listitems), 'ListSize', [dlgwidth, 300]);
end
if tf ~= 1 || isempty(indx), return; end

[~, methodid] = gui.i_pickscoremethod([], parentfig);
if isempty(methodid), return; end

sel = rows(indx);
n = numel(sel);
y = cell(n, 1); ttxt = cell(n, 1); posgc = cell(n, 1);
valid = false(n, 1);

fw = gui.myWaitbar(parentfig);
for k = 1:n
    % Brace-index to a string scalar rather than a 1-by-1 cell, so every
    % label reaching i_scoreNames and the plot has the same shape.
    ttxt{k} = string(T.ScoreType{sel(k)});
    % Skip a signature with too few expressed genes (which errors in
    % e_cellscores) instead of aborting the run.
    try
        [y{k}, ~, posgc{k}] = pkg.e_cellscores(sce.X, ...
            sce.g, sel(k), methodid, false);
    catch ME
        warning('Skipping score "%s": %s', ...
            string(T.ScoreType(sel(k))), ME.message);
        y{k} = [];
    end
    valid(k) = ~isempty(y{k});
end
gui.myWaitbar(parentfig, fw);

if ~any(valid)
    gui.myWarndlg(parentfig, ['No scores could be computed. The ' ...
        'selected signatures have too few expressed genes in this ' ...
        'dataset.']);
    y = {}; ttxt = {}; methodid = [];
    return;
end
if ~all(valid)
    gui.myWarndlg(parentfig, sprintf( ...
        '%d of %d score(s) skipped (too few expressed genes): %s', ...
        sum(~valid), n, ...
        strjoin(string(T.ScoreType(sel(~valid))), ', ')));
    y = y(valid); ttxt = ttxt(valid); posgc = posgc(valid);
end
posg = posgc{end};
end


function figfun = i_pickScatterPlotType(parentfig)
%I_PICKSCATTERPLOTTYPE Let the user choose heat or stem scatter.
%   Mirrors the chooser in scgeatoolApp/PushTool_ShowCellStatesClicked.
%   i_heatscatterfig and i_stemscatterfig take the same arguments, so the
%   choice is returned as a function handle. Empty means the user cancelled.

figfun = [];
switch gui.myQuestdlg(parentfig, 'Plot scatter plot type:', '', ...
        {'Heat', 'Stem'}, 'Heat')
    case 'Heat'
        figfun = @gui.i_heatscatterfig;
    case 'Stem'
        figfun = @gui.i_stemscatterfig;
end
end


function plotkind = i_pickComparisonPlotType(parentfig)
%I_PICKCOMPARISONPLOTTYPE Violin or bar for the group-comparison view.
%   Returns the plotkind string sc_uitabgrpfig_vioplot expects, or empty if
%   the user cancelled. Either way the figure keeps its "Switch to Bar Plot"
%   toolbar button, so the choice is not final.

plotkind = '';
switch gui.myQuestdlg(parentfig, 'Plot type:', '', ...
        {'Violin', 'Bar'}, 'Violin')
    case 'Violin'
        plotkind = 'violin';
    case 'Bar'
        plotkind = 'bar';
end
end


function choice = i_pickMultiScorePlotType(parentfig)
%I_PICKMULTISCOREPLOTTYPE Plot type for several scores at once.
%   'Heat' and 'Stem' give one tab per score via sc_uitabgrpfig_feaplot
%   (methodid 1 and 2 respectively), matching the single-score chooser.
%   'Score heatmap' is the scores-by-cells view this branch used to produce
%   unconditionally; it is still offered, but Heat is the default so both
%   choosers behave the same way.

choice = gui.myQuestdlg(parentfig, 'Plot type:', '', ...
    {'Heat', 'Stem', 'Score heatmap'}, 'Heat');
end


function [sce, saved] = i_saveScoresAsAttributes(sce, keepidx, ...
    y, ttxt, parentfig)
%I_SAVESCORESASATTRIBUTES Ask the user, then append score(s) to the SCE.
%
%   Handles both a single score (numeric vector in y) and several
%   (cell array of vectors). Names that collide with existing attributes are
%   made unique rather than silently dropped.
%
%   Scores computed on a subset of cells are not offered for saving at all -
%   see the keepidx guard below.

saved = false;

% Comparing two selected groups subsets the SCE, so the scores cover only
% those cells. A partial attribute is not worth storing, and pushing the
% subset SCE back to the app would discard every other cell, so skip the
% save entirely rather than asking a question with no good answer.
if ~all(keepidx), return; end

if iscell(y)
    values = y(:);
else
    values = {y};
end
values = values(~cellfun(@isempty, values));
if isempty(values), return; end

names = i_scoreNames(ttxt, numel(values));
if numel(names) ~= numel(values), return; end

if isscalar(values)
    prompt = sprintf('Save score "%s" to the cell attribute list?', names(1));
else
    prompt = sprintf('Save %d computed scores to the cell attribute list?', ...
        numel(values));
end

answer = gui.myQuestdlg(parentfig, prompt, '', ...
    {'Yes, save', 'No, skip'}, 'Yes, save');
if ~strcmp(answer, 'Yes, save'), return; end

for k = 1:numel(values)
    v = double(values{k}(:));
    if numel(v) ~= sce.NumCells
        warning('Skipping score "%s": %d values for %d cells.', ...
            names(k), numel(v), sce.NumCells);
        continue;
    end
    existing = string(sce.list_cell_attributes(1:2:end));
    thisname = names(k);
    if any(strcmp(thisname, existing))
        unique_names = matlab.lang.makeUniqueStrings([existing, thisname]);
        fprintf('Duplicate found: %s is renamed to %s.\n', ...
            thisname, unique_names(end));
        thisname = unique_names(end);
    end
    % Store the name as char. gui.i_select1state concatenates these into a
    % cell of char vectors, and a string scalar makes that cell mixed-type,
    % which uilistbox rejects ("'Items' must be a 1-by-N cell array of
    % character vectors or a string array").
    sce.list_cell_attributes = [sce.list_cell_attributes, {char(thisname), v}];
    saved = true;
end
end


function names = i_scoreNames(ttxt, n)
%I_SCORENAMES Normalise the assorted label shapes into an n-element string.
%   ttxt arrives as a char, a string scalar, a string array, or a cell of
%   any of those, depending on which score source produced it.

names = strings(0, 1);
if isempty(ttxt), return; end

if iscell(ttxt)
    names = arrayfun(@(k) string(ttxt{k}), (1:numel(ttxt))');
elseif ischar(ttxt)
    names = string(ttxt);   % a char row vector is one name, not one per char
else
    names = string(ttxt(:));
end
names = names(:);

% One label covering several scores: number them so they stay distinct.
if isscalar(names) && n > 1
    names = names + "_" + string((1:n)');
end
end

function txt = i_methodlabel(methodid)
% Human-readable name of the scoring algorithm, for the figure subtitle.
% gui.i_pickscoremethod returns the name outright for a known id and only
% opens a dialog for 0 or empty - so guard, or asking what the method was
% would ask the user to choose one all over again.
txt = '';
if isempty(methodid) || ~isscalar(methodid) || ~ismember(methodid, [1 2 3])
    return;
end
% The PMID tag earns its place in the chooser, not on an axis label.
txt = gui.i_stripcitation(gui.i_pickscoremethod(methodid));
end

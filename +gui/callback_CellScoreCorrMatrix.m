function callback_CellScoreCorrMatrix(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    [done, scorefile] = in_getscorfile;
    if ~done, return; end

    goodfile = false;
    if isfile(scorefile)
        try
            T = readtable(scorefile, 'Sheet', 'Sheet1', ...
                'ReadVariableNames', true);
            goodfile = true;
        catch ME
            disp(ME.message);
        end
    end
    if ~goodfile, return; end
    gsets = T.PositiveMarkers;
    
    answer = gui.myQuestdlg(FigureHandle, 'All cells, or select cells?', ...
        'Apply Analysis to All Cells?',{'All Cells','Select Cells'});
    switch answer
        case 'All Cells'
        case 'Select Cells'
            [thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle);
            if isempty(thisc), return; end
            [~, cL] = findgroups(string(thisc));
            [idx] = gui.i_selmultidlg(cL, [], FigureHandle);
            if isempty(idx), return; end

            sce = sce.selectcells(ismember(thisc, cL(idx)));
        otherwise
            return;
    end

    %{
    answer = gui.myQuestdlg(FigureHandle, 'In the heatmap, rearrange the order of gene programs to form clusters?', ...
        '');
    switch answer
        case 'Yes'
            docluster = true;
        case 'No'
            docluster = false;
        otherwise
            return;
    end
    %}


    fw = gui.myWaitbar(FigureHandle);
    [M, C] = pkg.e_cellscorecorrmat(sce.X, sce.g, gsets, 2, FigureHandle);
    labels = strrep(string(T.ScoreType),'_','\_');

    %if docluster
        t = clusterdata(C', maxclust=5);
        [~,idx]=sort(t);
        M2=M(idx,idx);
        labels2=labels(idx);
    %end
    % assignin("base","C",C);
    gui.myWaitbar(FigureHandle,fw);


    figure;
    isupper = logical(triu(ones(size(M)),0));
    M(isupper) = NaN;
    h = heatmap(M,'MissingDataColor','w');
    h.YDisplayLabels = labels;
    colormap(h, 'parula');
    title(h,'Gene programs in original order');
    %h.Colormap = flipud(h.Colormap);
    % imagesc(M);


    figure
    isupper = logical(triu(ones(size(M2)),0));
    M2(isupper) = NaN;
    h2 = heatmap(M2,'MissingDataColor','w');
    h2.YDisplayLabels = labels2;
    colormap(h2, 'parula');
    title(h2,'Gene programs are reordered to show any clusters');
    


    function [done, scorefile] = in_getscorfile
        done = false;
        pw1 = fileparts(mfilename('fullpath'));
        defaultscorefilename = 'cellscorecorrmat.xlsx';
        defaultscorefile = fullfile(pw1, '..', 'assets', 'CellScores', defaultscorefilename);
        preftagname =  'scorcorrmatfile';

        scorefile = getpref('scgeatoolbox', preftagname, defaultscorefile);
        if ~ispref('scgeatoolbox', preftagname)
            % if ~strcmp('Yes', gui.myQuestdlg(parentfig, 'Locate cellscorecorrmat.xlsx?')), return; end
            % [file, path] = uigetfile(defaultscorefilename, 'Select File');
            % if isequal(file, 0), return; end
            % scorefile = fullfile(path, file);
            % setpref('scgeatoolbox', preftagname, scorefile);
            % gui.myHelpdlg(parentfig, defaultscorefilename + " is located successfully.");
        else
            % scorefile = getpref('scgeatoolbox', preftagname);
            answer1 = gui.myQuestdlg(FigureHandle, sprintf('%s', scorefile), ...
                'Selected File', ...
                {'Use this', 'Use another', 'Cancel'}, 'Use this');
            if isempty(answer1), return; end
            switch answer1
                case 'Use this'
                    done = true;
                case 'Cancel'
                    return;
                case 'Use another'
                    % absolutePath = matlab.io.file.absolutePath(defaultscorefile);
                    % absolutePath=fileparts(which(defaultscorefile));
                    absolutePath = char(java.io.File(defaultscorefile).getCanonicalPath());
                    %javaFile = java.io.File(defaultscorefile);
                    %absolutePath = char(javaFile.getAbsolutePath());
                    [file, path] = uigetfile(defaultscorefilename, ...
                        'Select File', absolutePath);
                    if isequal(file, 0), return; end
                    scorefile = fullfile(path, file);
                    if isfile(scorefile)
                        setpref('scgeatoolbox', preftagname, scorefile);
                        gui.myHelpdlg(FigureHandle, defaultscorefilename + " is located successfully.");
                        done = true;
                    end
            end
        end

    end

end

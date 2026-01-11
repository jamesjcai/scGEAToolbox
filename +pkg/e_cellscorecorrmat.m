function [M, C] = e_cellscorecorrmat(X, g, gsets, methodid, ~)
    
    if nargin<4, methodid = 2; end
    if nargin<5, parentfig = []; end

    % preftagname = 'llapikeyenvfile';
    %{
    pw1 = fileparts(mfilename('fullpath'));    
    defaultscorefilename = 'cellscorecorrmat.xlsx';
    defaultscorefile = fullfile(pw1, '..', 'assets', 'CellScores', defaultscorefilename);
    goodfile = false;
    if isfile(defaultscorefile)
        try
            T = readtable(defaultscorefile, 'Sheet', 'Sheet1', ...
                'ReadVariableNames', true);
            goodfile = true;
        catch ME
            disp(ME.message);
        end
    end
    if ~goodfile, return; end
    gsets = T.PositiveMarkers;
   
% unique(strsplit(string(gsets{2}),','))

    preftagname = 'cellscorecorrmatfile';
    if ~ispref('scgeatoolbox', preftagname)
        if ~strcmp('Yes', gui.myQuestdlg(parentfig, 'Locate cellscores.xlsx?')), return; end
        [file, path] = uigetfile('cellscores.xlsx', 'Select File');
        if isequal(file, 0), return; end
        scorefile = fullfile(path, file);
        setpref('scgeatoolbox', preftagname, scorefile);
        gui.myHelpdlg(parentfig, "cellscores.xlsx is located successfully.");
    else
        scorefile = getpref('scgeatoolbox', preftagname);
    end
    %   [posg, ctselected] = gui.i_selectMSigDBGeneSets(speciestag, false, FigureHandle);
    %}

C = zeros(size(X, 2), length(gsets));

for k = 1:length(gsets)    
    tgsPos = unique(strsplit(string(gsets{k}),','));
    if methodid == 1
        [cs] = sc_cellscore_ucell(X, g, tgsPos);
    elseif methodid == 2
        [cs] = sc_cellscore_admdl(X, g, tgsPos);
    end
    C(:, k) = cs(:);
end
M = corr(C,'Type','Spearman');

% 
% 
% for k = 1:n
%     [y{k}, ~, posg] = pkg.e_cellscores(sce.X, sce.g, ...
%         indx2(k), methodid, false);
%     ttxt{k} = T.ScoreType(indx2(k));
% end
% 
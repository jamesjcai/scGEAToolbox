function [score, T, posg] = e_cellscores(X, genelist, ...
    scoretypeid, methodid, showwaitbar)
% Calcute predefined cell scores (marker list in cellscores.xlsx)
%
% see also: SC_CELLSCORE_UCELL, SC_CELLSCORE_ADMDL, SC_CELLCYCLESCORING

% scoretypeid

if nargin < 5, showwaitbar = true; end
if nargin < 4, methodid = []; end
if nargin < 3, scoretypeid = 0; end
if nargin < 2, genelist = []; end
if nargin < 1, X = []; end

score = [];

pw1 = fileparts(mfilename('fullpath'));
try
    cellscoresfile = fullfile(pw1, '..', 'assets', 'CellScores', ...
        'cellscores.xlsx');
    T = readtable(cellscoresfile, 'Sheet', 'Sheet1', ...
        'ReadVariableNames', true);
catch ME
    disp(ME.message);
    cellscoresfile = fullfile(pw1, '..', 'assets', ...
        'CellScores', 'cellscores.txt');
    T = readtable(cellscoresfile, 'Delimiter', '\t', ...
        'ReadVariableNames', true);
end

%T=sortrows(T,"ScoreType");

if ischar(scoretypeid) || isstring(scoretypeid)
    idx = find(matches(T.ScoreType, scoretypeid, 'IgnoreCase', true));
elseif isnumeric(scoretypeid)
    idx = scoretypeid;
end

if isempty(idx), return; end

idx = idx(1);    % in case there is duplicate in the score name list
if ~(idx <= height(T) && idx > 0 && idx == floor(idx))
    return;
end

scoretype = string(T.ScoreType(idx));

tgsPos = unique(pkg.i_str2genelist(T.PositiveMarkers(idx)));
tgsNeg = unique(pkg.i_str2genelist(T.NegativeMarkers(idx)));
if ~isempty(tgsNeg)
    tgsNeg = setdiff(tgsNeg, tgsPos);
end

% --------------------
posg = sort(tgsPos);
isexpressed = ismember(upper(posg), upper(genelist));
if sum(isexpressed)<2, error('Too few expressed genes (n < 2).'); end

answer = gui.i_pickscoremethod(methodid);
switch answer
    case 'AddModuleScore/Seurat'
        if showwaitbar, fw = gui.gui_waitbar; end
        try
            [score] = sc_cellscore_admdl(X, genelist, tgsPos, tgsNeg);
        catch ME
            if showwaitbar, gui.gui_waitbar(fw, true); end
            errordlg(ME.message);
            return;
        end
        if showwaitbar, gui.gui_waitbar(fw); end
    case 'UCell [PMID:34285779]'
        %[cs]=run.UCell(sce.X,sce.g,posg);

        if showwaitbar, fw = gui.gui_waitbar([], [], scoretype); end
        try
            [score] = sc_cellscore_ucell(X, genelist, tgsPos);
        catch ME
            if showwaitbar, gui.gui_waitbar(fw, true); end
            errordlg(ME.message);
            return;
        end
        if showwaitbar, gui.gui_waitbar(fw); end
    otherwise
        return;
end


fprintf('\n=============\n%s (%s)\n-------------\n', 'Genes', scoretype);
for k = 1:length(posg)
    if isexpressed(k)
        fprintf('%s*, ', posg(k));
    else
        fprintf('%s, ', posg(k));
    end
    if mod(k, 10) == 0 || k == length(posg)
        fprintf('\n');
    end
end
fprintf('=============\n*Expressed genes (n = %d)\n', ...
    sum(isexpressed));

% fprintf('\n=============\n%s\n-------------\n','Marker Genes');
% for k=1:length(posg)
%     fprintf('%s\n',posg(k));
% end
% fprintf('=============\n');

end

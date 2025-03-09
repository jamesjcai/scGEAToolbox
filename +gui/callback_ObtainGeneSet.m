function [posg, ttxt] = callback_ObtainGeneSet(~, ~)

FigureHandle = [];

selitems = {'Predefined Sets...', ...
    'MSigDB Signatures...', ...
    'PanglaoDB Cell Type Markers...', ...
    'TF Target Sets...'};

[indx1, tf1] = listdlg('PromptString', ...
    'Select a gene set from:', ...
    'SelectionMode', 'single', ...
    'ListString', selitems, ...
    'ListSize', [220, 300]);
if tf1 ~= 1, return; end

selecteditem = selitems{indx1};

switch selecteditem
    case 'Predefined Sets...'
        [~, T] = pkg.e_cellscores([], [], 0);
        listitems = T.ScoreType;
        [indx2, tf2] = listdlg('PromptString', 'Select Score', ...
            'SelectionMode', 'single', 'ListString', ...
            listitems, 'ListSize', [260, 300]);
        if tf2 ~= 1, return; end

        ttxt = string(T.ScoreType(indx2));
        posg = unique(strsplit(string(T.PositiveMarkers(indx2)), ...
            ','), 'stable');
        posg = posg(strlength(posg) > 0);

    case 'MSigDB Signatures...'
        stag = gui.i_selectspecies(2, true, FigureHandle);
        if isempty(stag), return; end
        try
            [posg, ttxt] = gui.i_selectMSigDBGeneSet(stag);
        catch ME
            errordlg(ME.message);
            return;
        end
        if isempty(posg) || isempty(ttxt), return; end
    case 'PanglaoDB Cell Type Markers...'
        stag = gui.i_selectspecies(2, true, FigureHandle);
        if isempty(stag), return; end

        oldpth = pwd;
        pw1 = fileparts(mfilename('fullpath'));
        pth = fullfile(pw1, '..', '+run', 'thirdparty', 'alona_panglaodb');
        cd(pth);

        markerfile = sprintf('marker_%s.mat', stag);
        if exist(markerfile, 'file')
            load(markerfile, 'Tm');
        else
            % Tw=readtable(sprintf('markerweight_%s.txt',stag));
            Tm = readtable(sprintf('markerlist_%s.txt', stag), ...
                'ReadVariableNames', false, 'Delimiter', '\t');
            % save(markerfile,'Tw','Tm');
        end
        cd(oldpth);

        ctlist = string(Tm.Var1);
        listitems = sort(ctlist);
        [indx, tf] = listdlg('PromptString', ...
            {'Select Class'}, ...
            'SelectionMode', 'single', ...
            'ListString', listitems, ...
            'ListSize', [220, 300]);
        if ~tf == 1, return; end
        ctselected = listitems(indx);
        % idx=find(matches(ctlist,ctselected));
        idx = matches(ctlist, ctselected);
        ctmarkers = Tm.Var2{idx};
        posg = string(strsplit(ctmarkers, ','));
        posg(strlength(posg) == 0) = [];
        ttxt = ctselected;

    case 'TF Target Sets...'
        [T] = pkg.e_gettflist;
        listitems = unique(T.tf);

        [indx2, tf2] = listdlg('PromptString', 'Select a transcription factor (TF)', ...
            'SelectionMode', 'single', 'ListString', ...
            listitems, 'ListSize', [220, 300]);
        if tf2 ~= 1, return; end

        ttxt = listitems{indx2};
        posg = string(T.target(string(T.tf) == ttxt));
    otherwise
        return;
end
% assignin('base','y',y);
% assignin('base','thisc',thisc);

if ~exist("posg", "var"), posg = []; end
if ~exist("ttxt", "var"), ttxt = []; end


posg = posg(:);
if nargout < 1
    inputdlg(ttxt, '', [15, 80], {char(posg(:))});
end

end

function sc_explorer(X, genelist, s, varargin)
%Single cell explorer
%
% pw1=fileparts(which(mfilename));
% pw1=fileparts(mfilename('fullpath'));
% pth=fullfile(pw1,'thirdparty/backup');
% addpath(pth);

p = inputParser;
addRequired(p, 'X', @isnumeric);
addRequired(p, 'genelist', @isstring);
addRequired(p, 's', @isnumeric);
addOptional(p, 'type', "unknown", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["unknown", "celltype", "marker", "pseudotime"]));
parse(p, X, genelist, s, varargin{:});
type = p.Results.type;

if strcmp(type, "unknown"), type = i_getoption; end

switch type
    case "celltype"
        MenuSelected1;
        return;
    case "marker"
        MenuSelected2;
        return;
    case "pseudotime"
        MenuSelected3;
        return;
    otherwise
        % return;
end

hFig = figure();
hAx = axes('Parent', hFig);
if size(s, 2) >= 5
    i_view5d(s);
    hFig.Position(3) = hFig.Position(3) * 2;
elseif size(s, 2) == 3
    scatter3(hAx, s(:, 1), s(:, 2), s(:, 3), 10);
elseif size(s, 2) == 2
    scatter(hAx, s(:, 1), s(:, 2), 10);
end

%tb = findall(hFig,'Type','uitoolbar');
tb = uitoolbar(hFig);
pt = uipushtool(tb, 'Separator', 'off');
[img, map] = imread(fileparts(mfilename('fullpath')), ...
    'resources', 'explorer1.gif');
ptImage = ind2rgb(img, map);
pt.CData = ptImage;
pt.Tooltip = 'Cell Type Explorer...';
pt.ClickedCallback = @MenuSelected1;

pt2 = uipushtool(tb, 'Separator', 'on');
[img, map] = imread(fileparts(mfilename('fullpath')), ...
    'resources', 'explorer2.gif');
ptImage = ind2rgb(img, map);
pt2.CData = ptImage;
pt2.Tooltip = 'Marker Gene Explorer...';
pt2.ClickedCallback = @MenuSelected2;

pt3 = uipushtool(tb, 'Separator', 'on');
[img, map] = imread(fileparts(mfilename('fullpath')), ...
    'resources', 'explorer3.gif');
ptImage = ind2rgb(img, map);
pt3.CData = ptImage;
pt3.Tooltip = 'Pseudotime Explorer...';
pt3.ClickedCallback = @MenuSelected3;

m = uimenu('Text', '&Explorers');
mitem1 = uimenu(m, 'Text', '&Cell Type Explorer...');
mitem2 = uimenu(m, 'Text', '&Marker Gene Explorer...');
mitem3 = uimenu(m, 'Text', '&Pseudotime Explorer...');
mitem1.Accelerator = 'C';
mitem2.Accelerator = 'M';
mitem3.Accelerator = 'P';
mitem1.MenuSelectedFcn = @MenuSelected1;
mitem2.MenuSelectedFcn = @MenuSelected2;
mitem3.MenuSelectedFcn = @MenuSelected3;
add_3dcamera;

%         cm = uicontextmenu(hFig);
%         m1 = uimenu(cm,'Text','Cell Type Explorer...');
%         m2 = uimenu(cm,'Text','Marker Gene Explorer...');
%         m3 = uimenu(cm,'Text','Pseudotime Explorer...');
%         % cm.ContextMenuOpeningFcn = @zzz;
%         hFig.ContextMenu = cm;
%         m1.MenuSelectedFcn = @MenuSelected1;
%         m2.MenuSelectedFcn = @MenuSelected2;
%         m3.MenuSelectedFcn = @MenuSelected3;


    function MenuSelected1(~, ~)
        answer = questdlg('Which species?', ...
            'Select Species', ...
            'Mouse', 'Human', 'Mouse');
        switch answer
            case 'Human'
                speciesx = "human";
            case 'Mouse'
                speciesx = "mouse";
            otherwise
                return;
        end
        answer = questdlg('Which algorithm?', ...
            'Select Method', ...
            'Alona', 'SingleR', 'Alona');
        switch answer
            case 'Alona'
                methodx = "alona";
            case 'SingleR'
                methodx = "singler";
            otherwise
                return;
        end
        gui.sc_celltypeexplorer(X, genelist, s, ...
            'species', speciesx, "method", methodx);
end

        function MenuSelected2(~, ~)
            gui.sc_markerexplorer(X, genelist, s);
    end
            function MenuSelected3(~, ~)
                gui.sc_pseudotimeexplorer(X, genelist, s);
        end

        end


            function [type] = i_getoption
            answer = questdlg('Select the type of explorer', ...
                'Select Explorer Type', ...
                'Cell type', 'Marker gene', 'Pseudotime', 'Cell type');
            switch answer
                case 'Cell type'
                    type = "celltype";
                case 'Marker gene'
                    type = "marker";
                case 'Pseudotime'
                    type = "pseudotime";
                otherwise
                    type = "unknown";
            end
        end
